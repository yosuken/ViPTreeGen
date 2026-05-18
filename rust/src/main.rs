//
// viptreegen-summary-pre / ViPTreeGen
//
//   Copyright: 2026 (C) Yosuke Nishimura
//   License: MIT
//
// Per-node parallel reader of per-(query)-node tblastx.out files. Replaces the
// legacy Ruby+SQL pipeline that spanned four pipeline steps:
//   01-3.cat_and_rename_split_tblastx
//   02-1.tblastx_filter
//   02-2.make_sbed_of_blast
//   02-3.make_summary_pre
// into a single binary that:
//   * lists per-node tblastx.out files (split-aware: <kind>/<qid>/blast/split/<i>/tblastx.out
//     when present, else <kind>/<qid>/blast/tblastx.out)
//   * applies qstart/qend shift = cutlen * (split_idx - 1)
//   * filters HSPs by pident > idt and alen >= alen_min
//   * converts to BED-like sub/que ranges (1-based inclusive)
//   * sweep-line interval merge per (lab, side), max identity & bitscore-per-bp per segment
//   * emits one summary_pre row per (kind, node, lab) ready for COPY into DuckDB
//
// CLI:
//   viptreegen-summary-pre --outdir DIR --kind {node|input} --cutlen N
//                          --idt F --alen N --threads N --output FILE
//
// Memory: O(threads * max_per_node_hsps * ~100B). Bench EVG (1811 seq, 12 threads,
// max 27k HSPs/node) peaks at ~280 MB. Scales to ~16 GB only beyond ~1M sequences.

use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Clone)]
struct Hit {
    split_key: String,
    sid: String,
    sub_s: i64,
    sub_e: i64,
    que_s: i64,
    que_e: i64,
    qlabel: String,
    identity: f64,
    bitscore: f64,
    evalue: f64,
}

struct Opts {
    outdir: PathBuf,
    kind: String,
    cutlen: i64,
    idt: f64,
    alen: i64,
    threads: usize,
    output: PathBuf,
}

fn parse_args() -> Opts {
    let mut opts = Opts {
        outdir: PathBuf::new(),
        kind: String::new(),
        cutlen: 100_000,
        idt: 30.0,
        alen: 30,
        threads: 0,
        output: PathBuf::new(),
    };
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--outdir"  => opts.outdir  = args.next().expect("--outdir").into(),
            "--kind"    => opts.kind    = args.next().expect("--kind"),
            "--cutlen"  => opts.cutlen  = args.next().expect("--cutlen").parse().unwrap(),
            "--idt"     => opts.idt     = args.next().expect("--idt").parse().unwrap(),
            "--alen"    => opts.alen    = args.next().expect("--alen").parse().unwrap(),
            "--threads" => opts.threads = args.next().expect("--threads").parse().unwrap(),
            "--output"  => opts.output  = args.next().expect("--output").into(),
            "-h" | "--help" => { print_help(); std::process::exit(0); }
            x => { eprintln!("unknown arg: {x}"); print_help(); std::process::exit(2); }
        }
    }
    if opts.outdir.as_os_str().is_empty() || opts.kind.is_empty() || opts.output.as_os_str().is_empty() {
        eprintln!("missing required args"); print_help(); std::process::exit(2);
    }
    if opts.kind != "node" && opts.kind != "input" {
        eprintln!("--kind must be node|input"); std::process::exit(2);
    }
    opts
}

fn print_help() {
    eprintln!("Usage: viptreegen-summary-pre --outdir DIR --kind {{node|input}} \\");
    eprintln!("                              [--cutlen N] [--idt F] [--alen N] \\");
    eprintln!("                              [--threads N] --output FILE");
}

/// Discover per-node tblastx.out files. Returns (split_idx, path) pairs.
/// split_idx is 1-based; for non-split nodes the single file is reported with split_idx=1.
fn collect_split_files(blast_dir: &Path) -> Vec<(i64, PathBuf)> {
    let split_dir = blast_dir.join("split");
    if split_dir.is_dir() {
        let mut v: Vec<(i64, PathBuf)> = Vec::new();
        if let Ok(rd) = fs::read_dir(&split_dir) {
            for e in rd.flatten() {
                let p = e.path();
                if !p.is_dir() { continue; }
                let idx: i64 = p.file_name()
                    .and_then(|n| n.to_str())
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                if idx == 0 { continue; }
                let f = p.join("tblastx.out");
                if f.is_file() { v.push((idx, f)); }
            }
        }
        v.sort_by_key(|(i, _)| *i);
        v
    } else {
        let f = blast_dir.join("tblastx.out");
        if f.is_file() { vec![(1, f)] } else { vec![] }
    }
}

fn process_node(
    qid: &str,
    files: &[(i64, PathBuf)],
    is_split: bool,
    seq_len: &HashMap<String, i64>,
    cutlen: i64,
    idt_min: f64,
    alen_min: i64,
    kind: &str,
) -> Vec<String> {
    // 1) load + filter + bed convert; apply q.start/q.end shift for split files
    let mut hits: Vec<Hit> = Vec::new();
    for (idx, path) in files {
        let shift: i64 = if is_split { cutlen * (idx - 1) } else { 0 };
        let file = match File::open(path) { Ok(f) => f, Err(_) => continue };
        for line in BufReader::new(file).lines() {
            let l = match line { Ok(s) => s, Err(_) => continue };
            if l.is_empty() { continue; }
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() < 12 { continue; }
            let pident: f64 = c[2].parse().unwrap_or(0.0);
            let alen:   i64 = c[3].parse().unwrap_or(0);
            if !(pident > idt_min && alen >= alen_min) { continue; }
            let qstart: i64 = c[6].parse::<i64>().unwrap_or(0) + shift;
            let qend:   i64 = c[7].parse::<i64>().unwrap_or(0) + shift;
            let sstart: i64 = c[8].parse().unwrap_or(0);
            let send:   i64 = c[9].parse().unwrap_or(0);
            let evalue: f64 = c[10].parse().unwrap_or(0.0);
            let bitscore: f64 = c[11].parse().unwrap_or(0.0);
            let ss = sstart.min(send) - 1;
            let se = sstart.max(send);
            let qs = qstart.min(qend) - 1;
            let qe = qstart.max(qend);
            let identity = round_even(pident * 10.0);
            let qlabel = format!("{}:{}-{}", qid, qs, qe);
            let split_idx = qs / cutlen + 1;
            let split_key = split_idx.to_string();
            hits.push(Hit {
                split_key,
                sid: c[1].to_string(),
                sub_s: ss + 1, sub_e: se,
                que_s: qs + 1, que_e: qe,
                qlabel, identity, bitscore, evalue,
            });
        }
    }
    if hits.is_empty() { return Vec::new(); }

    // 2) legacy ASC sort so that HashMap.insert per-range keeps last-write
    hits.sort_by(|a, b| {
        a.split_key.cmp(&b.split_key)
            .then(a.sid.cmp(&b.sid))
            .then(a.sub_s.cmp(&b.sub_s))
            .then(a.sub_e.cmp(&b.sub_e))
            .then(a.qlabel.cmp(&b.qlabel))
            .then(b.bitscore.partial_cmp(&a.bitscore).unwrap_or(std::cmp::Ordering::Equal))
            .then(a.evalue.partial_cmp(&b.evalue).unwrap_or(std::cmp::Ordering::Equal))
    });

    // 3) group by sid (lab)
    let mut by_lab: HashMap<String, Vec<&Hit>> = HashMap::new();
    for h in &hits { by_lab.entry(h.sid.clone()).or_default().push(h); }

    // 4) per-lab sweep-line
    type Row = (f64, f64, i64, i64, i64, i64, f64, f64, f64, f64, String);
    let mut rows: Vec<Row> = Vec::with_capacity(by_lab.len());
    for (lab, lab_hits) in &by_lab {
        let (sub_len, sub_idt, sub_scr) = sweepline(lab_hits, |h| (h.sub_s, h.sub_e));
        let (que_len, que_idt, que_scr) = sweepline(lab_hits, |h| (h.que_s, h.que_e));
        if sub_len == 0 || que_len == 0 { continue; }
        let qlen = *seq_len.get(qid).unwrap_or(&1);
        let slen = *seq_len.get(lab).unwrap_or(&1);
        rows.push((
            que_scr, sub_scr, qlen, slen, que_len, sub_len,
            round1(que_idt / que_len as f64 / 10.0),
            round1(sub_idt / sub_len as f64 / 10.0),
            round1(que_len as f64 * 100.0 / qlen as f64),
            round1(sub_len as f64 * 100.0 / slen as f64),
            lab.clone(),
        ));
    }

    // 5) ranks (deterministic tie-break by lab ASC)
    let mut by_sub: Vec<usize> = (0..rows.len()).collect();
    by_sub.sort_by(|&a, &b| rows[b].1.partial_cmp(&rows[a].1).unwrap()
                            .then_with(|| rows[a].10.cmp(&rows[b].10)));
    let mut sub_rank = vec![0usize; rows.len()];
    for (rank, &idx) in by_sub.iter().enumerate() { sub_rank[idx] = rank + 1; }
    let mut by_que: Vec<usize> = (0..rows.len()).collect();
    by_que.sort_by(|&a, &b| rows[b].0.partial_cmp(&rows[a].0).unwrap()
                            .then_with(|| rows[a].10.cmp(&rows[b].10)));
    let mut que_rank = vec![0usize; rows.len()];
    for (rank, &idx) in by_que.iter().enumerate() { que_rank[idx] = rank + 1; }

    let mut out = Vec::with_capacity(rows.len());
    for (i, r) in rows.iter().enumerate() {
        let cols = [
            kind.to_string(),
            qid.to_string(),
            que_rank[i].to_string(),
            sub_rank[i].to_string(),
            qid.to_string(),
            r.10.clone(),
            (round_even(r.0) as i64).to_string(),
            (round_even(r.1) as i64).to_string(),
            r.2.to_string(),
            r.3.to_string(),
            r.4.to_string(),
            r.5.to_string(),
            fmt1(r.6),
            fmt1(r.7),
            fmt1(r.8),
            fmt1(r.9),
        ];
        out.push(cols.join("\t"));
    }
    out
}

fn sweepline<F: Fn(&Hit) -> (i64, i64)>(hits: &[&Hit], range_of: F) -> (i64, f64, f64) {
    let mut by_range: HashMap<(i64, i64), (f64, f64)> = HashMap::with_capacity(hits.len());
    for &h in hits {
        let (s, e) = range_of(h);
        if s > e { continue; }
        let len = (e - s + 1) as f64;
        by_range.insert((s, e), (h.identity, h.bitscore / len));
    }
    if by_range.is_empty() { return (0, 0.0, 0.0); }
    let ranges: Vec<((i64, i64), (f64, f64))> = by_range.into_iter().collect();

    let mut points: Vec<i64> = Vec::with_capacity(ranges.len() * 2);
    for &((s, e), _) in &ranges { points.push(s); points.push(e); }
    points.sort_unstable();
    points.dedup();

    let mut total_len: i64 = 0;
    let mut total_idt: f64 = 0.0;
    let mut total_scr: f64 = 0.0;

    for i in 0..points.len() {
        let p = points[i];
        let mut max_idt = f64::NEG_INFINITY;
        let mut max_scr = f64::NEG_INFINITY;
        let mut found = false;
        for &((s, e), (idt, scr)) in &ranges {
            if s <= p && e >= p {
                if idt > max_idt { max_idt = idt; }
                if scr > max_scr { max_scr = scr; }
                found = true;
            }
        }
        if found {
            total_len += 1;
            total_idt += max_idt;
            total_scr += max_scr;
        }
        if i + 1 < points.len() {
            let p_b = points[i + 1];
            let glen = p_b - p - 1;
            if glen > 0 {
                let mut g_idt = f64::NEG_INFINITY;
                let mut g_scr = f64::NEG_INFINITY;
                let mut gfound = false;
                for &((s, e), (idt, scr)) in &ranges {
                    if s <= p && e >= p_b {
                        if idt > g_idt { g_idt = idt; }
                        if scr > g_scr { g_scr = scr; }
                        gfound = true;
                    }
                }
                if gfound {
                    total_len += glen;
                    total_idt += g_idt * glen as f64;
                    total_scr += g_scr * glen as f64;
                }
            }
        }
    }
    (total_len, total_idt, total_scr)
}

fn round1(x: f64) -> f64 { format!("{:.1}", x).parse().unwrap_or(x) }
fn fmt1(x: f64) -> String { format!("{:.1}", x) }

/// Round to nearest integer, half-to-even (banker's rounding, IEEE 754 default).
/// `f64::round()` uses half-away-from-zero; this matches Ruby's `Float#round`
/// (>=2.4) and `printf("%.0f", ...)` semantics so the pipeline is uniform.
fn round_even(x: f64) -> f64 {
    let floor = x.floor();
    let diff  = x - floor;
    if diff < 0.5 {
        floor
    } else if diff > 0.5 {
        floor + 1.0
    } else {
        // exact half: round to nearest even
        if (floor as i64) % 2 == 0 { floor } else { floor + 1.0 }
    }
}

fn load_lengths(all_len_path: &Path) -> HashMap<String, i64> {
    let mut m = HashMap::new();
    let f = File::open(all_len_path).expect("open all.len");
    for line in BufReader::new(f).lines().flatten() {
        let mut it = line.split('\t');
        let id  = it.next().unwrap_or("").to_string();
        let len = it.next().unwrap_or("0").parse().unwrap_or(0);
        if !id.is_empty() { m.insert(id, len); }
    }
    m
}

fn main() {
    let opts = parse_args();
    if opts.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(opts.threads).build_global().unwrap();
    }

    let t0 = Instant::now();
    let all_len = opts.outdir.join("cat/all/all.len");
    let seq_len = load_lengths(&all_len);
    eprintln!("[{}] {} sequences loaded from {} in {:.2}s",
              opts.kind, seq_len.len(), all_len.display(), t0.elapsed().as_secs_f64());

    let t1 = Instant::now();
    let kind_dir = opts.outdir.join(&opts.kind);
    let mut nodes: Vec<(String, bool, Vec<(i64, PathBuf)>)> = Vec::new();
    if let Ok(rd) = fs::read_dir(&kind_dir) {
        for e in rd.flatten() {
            let p = e.path();
            if !p.is_dir() { continue; }
            let qid = p.file_name().unwrap().to_string_lossy().to_string();
            let blast_dir = p.join("blast");
            if !blast_dir.is_dir() { continue; }
            let is_split = blast_dir.join("split").is_dir();
            let files = collect_split_files(&blast_dir);
            if !files.is_empty() { nodes.push((qid, is_split, files)); }
        }
    }
    nodes.sort_by(|a, b| a.0.cmp(&b.0));
    eprintln!("[{}] discovered {} nodes in {:.2}s",
              opts.kind, nodes.len(), t1.elapsed().as_secs_f64());

    let t2 = Instant::now();
    let outputs: Vec<Vec<String>> = nodes
        .par_iter()
        .map(|(qid, is_split, files)| {
            process_node(qid, files, *is_split, &seq_len, opts.cutlen, opts.idt, opts.alen, &opts.kind)
        })
        .collect();
    eprintln!("[{}] processed {} nodes in {:.2}s",
              opts.kind, outputs.len(), t2.elapsed().as_secs_f64());

    let t3 = Instant::now();
    let f = File::create(&opts.output).expect("create output");
    let mut w = BufWriter::new(f);
    let mut nrows: usize = 0;
    for lines in outputs {
        for line in lines {
            writeln!(w, "{}", line).expect("write");
            nrows += 1;
        }
    }
    w.flush().expect("flush");
    eprintln!("[{}] wrote {} summary_pre rows to {} in {:.2}s",
              opts.kind, nrows, opts.output.display(), t3.elapsed().as_secs_f64());
}
