"""The per-replicate deep-history path, from launcher to snakemake --config.

`cattle_baseline_seed` decides which ancestry a cattle run resumes from, and it
is the ONLY thing that distinguishes the four deep histories. If it fails to
reach the command line the run does not error -- it falls back to the config
YAML's 20250303 and the four histories collapse back into the one shared history
this work exists to remove, with nothing downstream reporting it.

Source-parsing rather than execution, in the style of `test_launcher_tables.py`
and `test_cattle_deep_tiers.py`: the wiring is what matters and it cannot be
exercised without a cluster.
"""

import os
import re

import pytest


SIM_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _read(name):
    with open(os.path.join(SIM_DIR, name)) as fh:
        return fh.read()


LAUNCHER = _read("submit_publication.sh")
CONTROLLER = _read("controller_publication.sbatch")
COMMON = _read("lib/publication_common.sh")


# ------------------------------------------------------- launcher -> sbatch

def test_the_launcher_resolves_the_deep_history_per_replicate():
    """Inside the replicate loop, not the category loop: E, L, F and G at one
    replicate share a history, and replicates do not."""
    assert "cattle_baseline_seed_of \"$R\"" in LAUNCHER


def test_only_cattle_gets_a_deep_history():
    """A config key the human path ignores is a key that can be wrong forever."""
    m = re.search(r'if \[\[ "\$SP" == "cattle" \]\]; then\s*\n\s*CB=\$\(cattle_baseline_seed_of',
                  LAUNCHER)
    assert m, "the deep-history lookup is not gated on species"


def test_the_seed_is_exported_to_the_controller():
    assert 'CB_SEED="${RUN_CB[$i]}"' in LAUNCHER
    assert 'SPECIES="${RUN_SPECIES[$i]}"' in LAUNCHER


def test_the_workdir_is_created_before_sbatch_writes_its_logs():
    """sbatch --output cannot create the directory, and the controller's own
    mkdir runs after slurmstepd has already tried to open the file. Harmless
    while every workdir survived a previous wave; not harmless for the 360
    replicate-6-and-up runs, whose directories have never existed."""
    i = LAUNCHER.index('mkdir -p "${RUN_WD[$i]}"')
    j = LAUNCHER.index('--output="${RUN_WD[$i]}/controller_%j.out"')
    assert i < j, "mkdir must precede the sbatch that writes into it"


def test_the_replicate_default_comes_from_the_table():
    """A literal here is how the manifest and the wave drift apart."""
    assert 'REPS="${REPS:-$(replicate_list)}"' in LAUNCHER
    assert "REPS=\"${REPS:-1 2 3 4 5}\"" not in LAUNCHER


# --------------------------------------------------------- the pre-flight

def test_the_preflight_iterates_the_wave_not_the_manifest():
    """The old version walked manifest rows and filtered them to the wave, so a
    run with NO row was never checked -- while printing a pass that counted the
    whole wave. With a 5-replicate manifest and a 20-replicate wave it reported
    120 runs healthy having verified 30."""
    assert 'for i in "${!RUN_DIR_NAMES[@]}"' in LAUNCHER
    assert "NO MANIFEST ROW" in LAUNCHER


def test_the_preflight_checks_the_marks_file():
    """rules/stage1_cattle_sel.smk requires the marks file AND the .full.ts. A
    present full file with a missing marks file falls through to a fresh SLiM
    run inside a 4 GB controller."""
    assert "m_marks" in LAUNCHER
    assert "stage1_marks_file" in LAUNCHER


def test_the_preflight_compares_the_deep_history_to_the_manifest():
    assert "DEEP HISTORY MISMATCH" in LAUNCHER


def test_the_preflight_reports_what_it_actually_checked():
    """The old pass line counted ${#RUN_IDS[@]} regardless of how many rows it
    had inspected, which is what made the hole invisible."""
    assert "checked $checked of ${#RUN_IDS[@]}" in LAUNCHER


# ------------------------------------------------- controller -> snakemake

def test_the_controller_refuses_a_cattle_run_without_a_deep_history():
    assert "a cattle run needs CB_SEED" in CONTROLLER
    assert "exit 1" in CONTROLLER.split("a cattle run needs CB_SEED")[1][:600]


def test_the_seed_reaches_both_snakemake_invocations():
    """The --unlock call and the real call. A stale claim takes the unlock path
    first, and a mismatch between the two would unlock one DAG and run another."""
    calls = [l for l in CONTROLLER.splitlines() if '"$STAGE1_CFG"' in l]
    assert len(calls) == 2, calls
    for line in calls:
        assert "$CB_CFG" in line


def test_the_controller_forbids_building_a_deep_history():
    """cattle_baseline_search_dirs=[] , not the deep-history directory. Stage 1
    is always adopted here; with the directory listed, a missing stage-1 file
    would start a ~2 h SLiM run inside a controller sized for symlinks."""
    assert "cattle_baseline_search_dirs=[]" in CONTROLLER


def test_the_deep_history_config_is_a_single_space_separated_token_list():
    """EXTRA_CONFIG and CB_CFG are expanded UNQUOTED so each key=value becomes
    its own --config argument. A token containing a space would silently become
    two arguments."""
    vals = re.findall(r'CB_CFG="([^"]*)"', CONTROLLER)
    assert vals, "no CB_CFG assignment found"
    body = [v for v in vals if v]          # the empty init is the other one
    assert len(body) == 1, vals
    for tok in body[0].split():
        assert "=" in tok and not tok.endswith("=")


# ------------------------------------------------------------ the library

def test_the_deep_histories_live_outside_stage1_inputs():
    """One directory holding both would put four ep7 handoffs next to the
    stage-1 files, and search_dirs._seed_collision_check cannot catch a wrong
    pick among them: its glob rewrites seed_<N> to sd*, which matches no
    seed_-named file."""
    assert 'DEEP_HISTORIES="${DEEP_HISTORIES:-$PUBLISH_ROOT/deep_histories}"' in COMMON


def test_the_library_stays_bash_3_2_clean():
    """The test suite sources this file, and macOS ships bash 3.2.57 -- no
    associative arrays, no ${x,,}."""
    assert "declare -A" not in COMMON
    assert not re.search(r"\$\{[A-Za-z_][A-Za-z0-9_]*,,\}", COMMON)


# ------------------------------------------- the deep-history builder itself

REBUILD = _read("rebuild_cattle_deep_history.sh")


def test_the_seed_reaches_snakemake():
    """THE BUG. The seed was read out of the configfile and used only to spell
    filenames for --rescale/--check; the --config list carried workdir and
    publishdir and nothing else. So the script could build exactly one deep
    history -- the one its configfile named -- and asking it for another would
    have silently rebuilt 20250303 under a different name."""
    wrap = REBUILD[REBUILD.index("--wrap="):]
    assert wrap.count("stage1_seed=$CB_SEED") == 2, \
        "both the --unlock and the real invocation must carry the seed"


def test_the_seed_is_accepted_from_the_environment():
    assert 'CB_SEED="${CB_SEED:-$(awk' in REBUILD


def test_the_workdir_is_seed_namespaced():
    """Filenames carry the seed; the workdir did not, and everything else under
    it is shared -- .snakemake locks and iocache, slurm_logs, params/, and the
    three rules' workdir-relative log: paths. Concurrent builds would collide."""
    m = re.search(r'^WD="\$\{WD:-([^"]*)\}"', REBUILD, re.M)
    assert m, "no WD default found"
    assert "$CB_SEED" in m.group(1), m.group(1)


def test_concurrent_builds_are_distinguishable_in_the_queue():
    assert '--job-name="deep_history_$CB_SEED"' in REBUILD


# ------------------------------------------------ the stage-1 regeneration

REGEN = _read("regenerate_stage1.sh")


def test_regen_builds_every_replicate_the_table_defines():
    """A literal 20-ID list would silently skip replicates 6-20."""
    assert "IDS=\"${IDS:-$(_default_ids)}\"" in REGEN
    assert "for r in $(replicate_list)" in REGEN


def test_regen_passes_the_per_replicate_deep_history():
    """Without the explicit seed a cattle replicate takes the config YAML's
    20250303 while its filename claims _cb{other} -- a stage-1 file that lies
    about its own ancestry, which nothing downstream could detect."""
    assert 'cattle_baseline_seed=$CB' in REGEN
    assert 'cattle_baseline_seed_of "${ID:1}"' in REGEN


def test_regen_looks_for_handoffs_in_the_deep_history_directory():
    assert "cattle_baseline_search_dirs=['$DEEP_HISTORIES']" in REGEN
    assert "seed_20250303.ep7.ts" not in REGEN, "a hardcoded handoff is left"


def test_regen_checks_every_handoff_the_wave_needs():
    assert "MISSING_HANDOFF" in REGEN


def test_regen_staleness_compares_against_the_right_handoff():
    """One reference would compare every replicate against block 1's handoff and
    report blocks 2-4 stale."""
    assert 'HANDOFF_REF="$DEEP_HISTORIES/$(deep_history_handoff_of "${ID:1}")"' in REGEN


# --------------------------------------------- relaunching a truncated wave

def test_a_partially_finished_run_is_resumed_not_refused():
    """`keep-going: True` lets a controller exit COMPLETED having lost one panel
    to a walltime kill, leaving 1 of 2 .enloc.sig.out. That run is exactly the
    one that must be relaunched, and the old predicate ("holds ANY stage-4
    output -> refuse the whole wave") refused it -- which made a truncated wave
    unfinishable without hand-editing the launcher. Block 2 came back with 30 of
    120 runs in that state."""
    assert "RESUME" in LAUNCHER
    assert "ENLOC_PER_RUN=2" in LAUNCHER
    assert "refusing to overwrite finished output" not in LAUNCHER


def test_completeness_covers_stage_five_not_just_stage_four():
    """A run can have both .enloc.sig.out and still be short a stage-5 panel --
    7 of block 2's 30 incomplete runs read `glmleads=2/want3` with no enloc
    complaint. Skipping on enloc alone would leave those unfinishable."""
    assert "GLM_PER_RUN=3" in LAUNCHER
    i = LAUNCHER.index("GLM_HITS=(")
    assert '"$g" -ge "$GLM_PER_RUN"' in LAUNCHER[i:i + 500]


def test_a_finished_run_is_skipped_not_refused():
    """Same reasoning as the live-lock skip: refusing the whole arm because part
    of it is already done is the wrong answer to a relaunch."""
    i = LAUNCHER.index("already complete (enloc")
    assert "skipped_done" in LAUNCHER[i - 400:i + 400]
    assert "FAIL=1" not in LAUNCHER[i - 300:i + 300]


def test_the_enloc_count_cannot_silently_exit_the_launcher():
    """`set -euo pipefail` + `ls | wc -l` = a silent exit when a run has ZERO
    enloc files: ls exits non-zero, pipefail propagates, the assignment fails
    and -e kills the script mid-loop having submitted nothing. That happened --
    the resume wave logged its SKIP/RESUME lines and then just stopped. Counted
    with a glob array instead; an unmatched glob stays literal so -f is false."""
    i = LAUNCHER.index("ENLOC_HITS=(")
    window = LAUNCHER[i - 200:i + 400]
    assert "ls -1" not in window
    assert "wc -l" not in window
    assert 'for f in "${ENLOC_HITS[@]}"' in window


# ------------------------------------- what counts as "already in progress"

def _keep_idx_loop():
    """The launcher's in-progress decision, isolated from the rest of the file.

    Scoped deliberately: `.snakemake/locks` may legitimately appear elsewhere
    (a diagnostic, a comment), and what this section asserts is what the
    DECISION reads.
    """
    i = LAUNCHER.index("declare -a KEEP_IDX=()")
    j = LAUNCHER.index("if (( skipped_locked ))", i)
    return LAUNCHER[i:j]


def test_in_progress_is_decided_by_the_claim_not_the_lock():
    """THE BUG that left 15 runs unfinishable.

    The loop globbed `.snakemake/locks/*` and called any hit "a controller
    already holds this workdir". But a snakemake killed at walltime dies without
    releasing its lock, while the controller's own `trap ... EXIT` still fires and
    removes the claim -- leaving a lock and no liveness record. Every relaunch
    then refused the run. Arm causal_maf001_paired reps 16-20 finished 15/30, its
    two passes submitted 1 job between them, and all 15 stragglers had locks=2,
    no claim, and no running job. The lock says snakemake once ran here; only the
    claim says whether anything runs here now.
    """
    loop = _keep_idx_loop()
    assert ".snakemake/locks" not in loop, \
        "the in-progress decision is reading the lock again"
    assert ".controller_claim/jobid" in loop


def test_only_a_running_or_pending_owner_blocks_a_relaunch():
    """A claim whose owner has finished, failed or been purged from slurmctld is
    not a live controller. Anything less specific than these two states makes a
    dead owner block its workdir permanently."""
    loop = _keep_idx_loop()
    assert 'scontrol show job "$OWNER"' in loop
    assert '"$STATE" == "JobState=RUNNING"' in loop
    assert '"$STATE" == "JobState=PENDING"' in loop


def test_the_launcher_and_the_controller_agree_on_liveness():
    """Two different answers to "who holds this workdir" is how the original
    double-controller incident happened. Same primitive, same states, both
    sides."""
    for src in (LAUNCHER, CONTROLLER):
        assert 'JobState=RUNNING' in src and 'JobState=PENDING' in src
        assert 'scontrol show job' in src


def test_the_skip_message_names_the_job_holding_the_workdir():
    """"a controller already holds this workdir" with no job id is what made this
    take a manual scratch crawl to diagnose."""
    loop = _keep_idx_loop()
    assert '"$OWNER"' in loop[loop.index("printf"):]


# --------------------------------------------- releasing an orphaned lock

def test_the_controller_unlocks_a_lock_left_under_a_fresh_claim():
    """The other half of the same bug. With the launcher fixed, a relaunch onto
    one of those workdirs takes a FRESH claim (no claim directory existed), so
    the old `CLAIMED == stale` guard skipped the unlock -- and snakemake exited
    on a stale lock it was entitled to release.

    Safe because a fresh claim is itself the proof: a live controller always
    holds its claim, so it would have sent this job out through REFUSING above.
    """
    line = [l for l in CONTROLLER.splitlines()
            if l.startswith('if [[ "$CLAIMED" == "stale" ]]')]
    assert len(line) == 1, line
    assert 'compgen -G "$WD/.snakemake/locks/*"' in line[0]


def test_the_unlock_is_still_conditional():
    """Unconditional unlocking is the original sin here: it cannot tell a lock
    held by a dead process from one held by a live one, and it released the live
    one. The condition must remain."""
    lines = CONTROLLER.splitlines()
    # the invocation, not the several comments that discuss it
    cmd = next(i for i, l in enumerate(lines)
               if "$SNAKEMAKE" in l and "--unlock" in l)
    guard = next(i for i in range(cmd - 1, -1, -1) if lines[i].startswith("if [["))
    assert '"$CLAIMED" == "stale"' in lines[guard]
