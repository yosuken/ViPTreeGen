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
#   - ruby (>= 2.0), gnu parallel
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
duckdb --version
ruby --version
test -x rust/target/release/viptreegen-summary-pre \
	|| fail "rust binary not built (run: cargo build --release --manifest-path rust/Cargo.toml)"

# ----- mode coverage ------------------------------------------------------------
for mode in tblastx mmseqs-tblastx blastn; do
	step "mode=$mode (normal, --notree)"
	out="$TEST_TMP/$mode"
	./ViPTreeGen --notree --ncpus "$NCPUS" --mode "$mode" testdata/ssDNA.prok.8.fasta "$out"
	test -f "$out/result/all.sim.matrix" || fail "$mode: result/all.sim.matrix missing"
done

# ----- numerical regression: golden matrix for --mode tblastx -------------------
# tblastx is deterministic across runs given the same BLAST+ version, so we ship
# a frozen matrix under testdata/expected/ and diff it bit-exact. If this fails
# after a BLAST+ version bump, regenerate with:
#   ./ViPTreeGen --notree --mode tblastx testdata/ssDNA.prok.8.fasta /tmp/golden \
#     && cp /tmp/golden/result/all.{sim,dist}.matrix testdata/expected/
# and document the BLAST+ version in the commit message.
step "golden matrix diff (tblastx all.sim.matrix bit-exact)"
diff -u testdata/expected/ssDNA.prok.8.tblastx.sim.matrix  "$TEST_TMP/tblastx/result/all.sim.matrix"  || fail "tblastx sim.matrix differs from golden"
diff -u testdata/expected/ssDNA.prok.8.tblastx.dist.matrix "$TEST_TMP/tblastx/result/all.dist.matrix" || fail "tblastx dist.matrix differs from golden"

# ----- label normalization ------------------------------------------------------
step "long-name normalization (exercises ParseFastaEntries)"
./ViPTreeGen --notree --ncpus "$NCPUS" testdata/ssDNA.prok.8.long_name.fasta "$TEST_TMP/longname"
test -f "$TEST_TMP/longname/result/all.sim.matrix" || fail "longname: result/all.sim.matrix missing"

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
