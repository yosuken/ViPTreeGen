#!/usr/bin/env bash
#
# ViPTreeGen smoke tests. Runs the same suite that CI runs (.github/workflows/test.yml)
# so a successful local run is a strong signal CI will pass too.
#
# Usage:
#   bash test.sh                  # default: 2 threads, all sub-tests
#   NCPUS=8 bash test.sh          # override thread count
#   SKIP_TREE=1 bash test.sh      # skip the R/ape/phangorn tree generation test
#   TEST_TMP=/some/path bash test.sh   # override scratch dir (default /tmp/vpt_test)
#
# Prerequisites (must be on PATH):
#   - tblastx, blastn, makeblastdb  (BLAST+; 'module load blast+' on JAMSTEC HPC)
#   - mmseqs                        ('module load mmseqs2')
#   - duckdb CLI
#   - ruby (>= 3.0), gnu parallel
#   - R + 'ape' + 'phangorn'        (only when SKIP_TREE is unset)
#   - rust/target/release/viptreegen-summary-pre  (build with `cargo build --release --manifest-path rust/Cargo.toml`)
#

set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

NCPUS="${NCPUS:-2}"
TEST_TMP="${TEST_TMP:-/tmp/vpt_test}"
SKIP_TREE="${SKIP_TREE:-}"

# Clean previous test artifacts (idempotent).
rm -rf "$TEST_TMP"
mkdir -p "$TEST_TMP"

step() {
	echo ""
	echo -e "\e[1;36m=== [test.sh] $* ===\e[0m"
}

fail() {
	echo -e "\e[1;31m[test.sh] FAILED: $*\e[0m" >&2
	exit 1
}

# ----- tool versions ------------------------------------------------------------
step "versions"
./ViPTreeGen --version
blastn -version | head -1
tblastx -version | head -1
mmseqs version | head -1
lastal --version | head -1
duckdb --version
ruby --version
test -x rust/target/release/viptreegen-summary-pre \
	|| fail "rust binary not built (run: cargo build --release --manifest-path rust/Cargo.toml)"

# ----- mode coverage ------------------------------------------------------------
for mode in tblastx mmseqs-tblastx blastn last; do
	step "mode=$mode (normal, --notree)"
	out="$TEST_TMP/$mode"
	./ViPTreeGen --notree --ncpus "$NCPUS" --mode "$mode" testdata/ssDNA.prok.8.fasta "$out"
	test -f "$out/result/all.sim.matrix" || fail "$mode: result/all.sim.matrix missing"
done

# ----- numerical regression: golden matrices for all four modes ----------------
# Each of tblastx / mmseqs-tblastx / blastn / last is deterministic given a fixed
# engine version + the same input, so we ship frozen matrices under
# testdata/expected/ and bit-exact diff against them. If this fails after a
# BLAST+ / MMseqs2 / LAST version bump, regenerate with:
#   for m in tblastx mmseqs-tblastx blastn last; do
#     ./ViPTreeGen --notree --mode "$m" testdata/ssDNA.prok.8.fasta "/tmp/g_$m" \
#       && cp /tmp/g_$m/result/all.{sim,dist}.matrix \
#          "testdata/expected/ssDNA.prok.8.$m.sim.matrix" \
#          "testdata/expected/ssDNA.prok.8.$m.dist.matrix" 2>/dev/null
#   done
# and document the new engine versions in the commit message.
for mode in tblastx mmseqs-tblastx blastn last; do
	step "golden matrix diff (--mode $mode bit-exact)"
	diff -u "testdata/expected/ssDNA.prok.8.$mode.sim.matrix"  "$TEST_TMP/$mode/result/all.sim.matrix"  || fail "$mode sim.matrix differs from golden"
	diff -u "testdata/expected/ssDNA.prok.8.$mode.dist.matrix" "$TEST_TMP/$mode/result/all.dist.matrix" || fail "$mode dist.matrix differs from golden"
done

# ----- label normalization ------------------------------------------------------
step "long-name normalization (exercises ParseFastaEntries)"
./ViPTreeGen --notree --ncpus "$NCPUS" testdata/ssDNA.prok.8.long_name.fasta "$TEST_TMP/longname"
test -f "$TEST_TMP/longname/result/all.sim.matrix" || fail "longname: result/all.sim.matrix missing"

# ----- gzip input ---------------------------------------------------------------
step "gzip-compressed input (auto-detected by magic bytes)"
gzip -c testdata/ssDNA.prok.8.fasta > "$TEST_TMP/ssDNA.prok.8.fasta.gz"
./ViPTreeGen --notree --ncpus "$NCPUS" --mode blastn "$TEST_TMP/ssDNA.prok.8.fasta.gz" "$TEST_TMP/gz"
# gz input must produce the exact same matrix as the plain golden
diff -u "testdata/expected/ssDNA.prok.8.blastn.sim.matrix" "$TEST_TMP/gz/result/all.sim.matrix" \
	|| fail "gzip input: sim.matrix differs from plain golden"

# ----- 2D mode ------------------------------------------------------------------
step "2D mode (query x input)"
./ViPTreeGen --notree --2D testdata/1.fasta --ncpus "$NCPUS" testdata/2.fasta "$TEST_TMP/2D"
test -f "$TEST_TMP/2D/result/2D.sim.matrix" || fail "2D: result/2D.sim.matrix missing"
test -f "$TEST_TMP/2D/result/top10.sim.list" || fail "2D: result/top10.sim.list missing"

# ----- with-reference mode (--ref-duckdb) ---------------------------------------
step "with-reference mode (build ref + reuse)"
# split testdata/1.fasta into ref (first 10) + input (last 5)
awk 'BEGIN{n=0} /^>/{n++} { if (n<=10) print > "'"$TEST_TMP"'/ref.fasta"; else print > "'"$TEST_TMP"'/input.fasta" }' testdata/1.fasta
./ViPTreeGen --notree --ncpus "$NCPUS" "$TEST_TMP/ref.fasta" "$TEST_TMP/ref_run"
./ViPTreeGen --notree --ncpus "$NCPUS" --ref-duckdb "$TEST_TMP/ref_run/run.duckdb" "$TEST_TMP/input.fasta" "$TEST_TMP/ref_use"
test -f "$TEST_TMP/ref_use/result/all.sim.matrix" || fail "ref mode: result/all.sim.matrix missing"
# matrix must be 15x15 (10 ref + 5 input) = 16 lines (header + 15 ids)
rows=$(wc -l < "$TEST_TMP/ref_use/result/all.sim.matrix")
[ "$rows" -eq 16 ] || fail "ref mode: expected 16 rows in all.sim.matrix, got $rows"

# ----- --mmseqs-tblastx-chunk-size: chunking is deterministic ------------------
# Force the 8-seq input through 3 query chunks and verify the resulting matrix
# matches the 1-chunk golden bit-exact. This protects the chunking machinery
# (ChunkFasta, per-chunk step_done, cat -> m8_split) from drift.
step "mmseqs-tblastx chunking: chunk_size=3 (-> 3 chunks) matches default golden"
CHUNK_DIR="$TEST_TMP/chunk3"
rm -rf "$CHUNK_DIR"
./ViPTreeGen --notree --ncpus "$NCPUS" --mmseqs-tblastx-chunk-size 3 \
             testdata/ssDNA.prok.8.fasta "$CHUNK_DIR"
diff -u testdata/expected/ssDNA.prok.8.mmseqs-tblastx.sim.matrix \
        "$CHUNK_DIR/result/all.sim.matrix" \
	|| fail "chunk_size=3 result differs from default (1-chunk) golden"
# Sanity: we should have step_done markers for 3 chunks
chunks=$(duckdb -noheader -list "$CHUNK_DIR/run.duckdb" \
	"SELECT COUNT(*) FROM run_metadata WHERE key LIKE 'step_done:01-2.chunk:%'")
[ "$chunks" = "3" ] || fail "expected 3 chunk step_done markers, got $chunks"

step "mmseqs-tblastx chunking: chunk_size=0 (disable) matches default golden"
NO_CHUNK_DIR="$TEST_TMP/nochunk"
rm -rf "$NO_CHUNK_DIR"
./ViPTreeGen --notree --ncpus "$NCPUS" --mmseqs-tblastx-chunk-size 0 \
             testdata/ssDNA.prok.8.fasta "$NO_CHUNK_DIR"
diff -u testdata/expected/ssDNA.prok.8.mmseqs-tblastx.sim.matrix \
        "$NO_CHUNK_DIR/result/all.sim.matrix" \
	|| fail "chunk_size=0 result differs from default golden"

# ----- --resume: simulated mid-flight kill --------------------------------------
# Simulate a kill mid-pipeline by deleting the last step's step_done marker, then
# re-run with --resume and verify the matrix still matches the golden afterwards.
step "--resume: re-run from a partial state, expect bit-exact final matrix"
RESUME_DIR="$TEST_TMP/resume_test"
rm -rf "$RESUME_DIR"
./ViPTreeGen --notree --ncpus "$NCPUS" --mode tblastx testdata/ssDNA.prok.8.fasta "$RESUME_DIR"
# Remove the final step's step_done marker AND the final matrix to force re-run of 03-1.
duckdb "$RESUME_DIR/run.duckdb" \
	"DELETE FROM run_metadata WHERE key = 'step_done:03-1.make_matrix'"
rm -f "$RESUME_DIR/result/all.sim.matrix" "$RESUME_DIR/result/all.dist.matrix"
./ViPTreeGen --resume --notree --ncpus "$NCPUS" --mode tblastx testdata/ssDNA.prok.8.fasta "$RESUME_DIR"
# After resume the final matrix should match the golden bit-exact.
diff -u testdata/expected/ssDNA.prok.8.tblastx.sim.matrix "$RESUME_DIR/result/all.sim.matrix" \
	|| fail "--resume final matrix differs from golden"

step "--resume: rejects parameter change vs prior run"
if ./ViPTreeGen --resume --notree --mode mmseqs-tblastx --ncpus "$NCPUS" \
                testdata/ssDNA.prok.8.fasta "$RESUME_DIR" 2>/dev/null; then
	fail "--resume should reject differing --mode (prev=tblastx, current=mmseqs-tblastx)"
fi

step "--resume: tblastx batch-level filter (SearchCmd tempfile pattern)"
# Verify SearchCmd writes the .tmp + mv pattern (not direct -out file).
batch_line=$(head -1 "$RESUME_DIR/batch/01-1"/*)
echo "$batch_line" | grep -q -- '-out .*\.tmp .*&& mv .*\.tmp' \
	|| fail "SearchCmd should emit '-out FILE.tmp ... && mv FILE.tmp FILE'; got: $batch_line"
# Simulate a real "01-2 killed mid-flight" state: 01-1 completed (its step_done is
# preserved), 01-2 onwards were killed (step_done cleared). Resume should re-enter
# 01-2, run the batch filter (no pre-existing tblastx.out -> 0 already done, N to
# run), and produce a matrix that still matches the tblastx golden bit-exact.
BATCH_DIR="$TEST_TMP/batch_resume"
rm -rf "$BATCH_DIR"
./ViPTreeGen --notree --ncpus "$NCPUS" --mode tblastx testdata/ssDNA.prok.8.fasta "$BATCH_DIR"
duckdb "$BATCH_DIR/run.duckdb" \
	"DELETE FROM run_metadata WHERE key LIKE 'step_done:01-2.%' \
	  OR key LIKE 'step_done:02-%' OR key LIKE 'step_done:03-%'"
log=$(./ViPTreeGen --resume --notree --ncpus "$NCPUS" --mode tblastx testdata/ssDNA.prok.8.fasta "$BATCH_DIR" 2>&1)
echo "$log" | grep -q "batch entries already done" \
	|| fail "expected resume to print '… batch entries already done, … to (re-)run' (log: $(echo \"$log\" | tail -5))"
diff -u testdata/expected/ssDNA.prok.8.tblastx.sim.matrix "$BATCH_DIR/result/all.sim.matrix" \
	|| fail "batch-resume matrix differs from golden"

step "--resume: rejects when output dir does not exist"
if ./ViPTreeGen --resume --notree --mode tblastx --ncpus "$NCPUS" \
                testdata/ssDNA.prok.8.fasta "$TEST_TMP/nonexistent_resume_dir" 2>/dev/null; then
	fail "--resume should reject when output dir does not exist"
fi

# ----- validation: invalid --mode rejected --------------------------------------
step "reject invalid --mode (mmseqs-blastn was removed; see doc/mmseqs-blastn.md)"
if ./ViPTreeGen --notree --mode mmseqs-blastn testdata/ssDNA.prok.8.fasta "$TEST_TMP/should_fail_a" 2>/dev/null; then
	fail "--mode mmseqs-blastn should have been rejected"
fi

# ----- validation: --2D + --ref-duckdb mutex ------------------------------------
step "reject conflicting --2D + --ref-duckdb"
touch "$TEST_TMP/dummy.duckdb"
if ./ViPTreeGen --notree --2D testdata/1.fasta --ref-duckdb "$TEST_TMP/dummy.duckdb" \
                testdata/2.fasta "$TEST_TMP/should_fail_b" 2>/dev/null; then
	fail "--2D + --ref-duckdb conflict should have been rejected"
fi

# ----- DuckDB schema sanity -----------------------------------------------------
step "run.duckdb schema sanity (schema_version=2.0, sequences.seq populated)"
duckdb "$TEST_TMP/mmseqs-tblastx/run.duckdb" '.tables'
duckdb "$TEST_TMP/mmseqs-tblastx/run.duckdb" \
	"SELECT key, value FROM run_metadata WHERE key IN ('schema_version','mode','search_mode')"
empty_seq=$(duckdb -noheader -list "$TEST_TMP/mmseqs-tblastx/run.duckdb" \
	"SELECT COUNT(*) FROM sequences WHERE seq IS NULL OR LENGTH(seq) = 0")
[ "$empty_seq" = "0" ] || fail "sequences.seq has $empty_seq empty rows (schema 2.0 contract broken)"

# ----- tree generation (slow, requires R + ape + phangorn) ----------------------
if [ -n "$SKIP_TREE" ]; then
	step "tree generation: SKIPPED (SKIP_TREE=$SKIP_TREE)"
else
	step "tree generation (exercises R + ape + phangorn)"
	./ViPTreeGen --ncpus "$NCPUS" --mode mmseqs-tblastx testdata/ssDNA.prok.8.fasta "$TEST_TMP/tree"
	test -f "$TEST_TMP/tree/result/all.sim.matrix" || fail "tree: result/all.sim.matrix missing"
	newicks=$(ls "$TEST_TMP/tree/result/"*.newick 2>/dev/null | wc -l)
	[ "$newicks" -gt 0 ] || fail "tree: no .newick file produced"
fi

# ----- done ---------------------------------------------------------------------
echo ""
echo -e "\e[1;32m=== [test.sh] all smoke tests passed ===\e[0m"
echo "Scratch outputs left under: $TEST_TMP"
