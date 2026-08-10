"""The cattle (Holstein-Friesian) epoch schedule, as backwards-time size changes.

WHY THIS FILE EXISTS

The demography is a pair of literal vectors that SLiM applies forward, one
`setSubpopulationSize()` per epoch:

    // farm_selection.slim:29-30, duplicated verbatim at
    // farm_selection_from_ep8.slim:171-172
    defineConstant("ep_sizes",   c(0, 17000, 10000, 7000, 3500, 2500, 2000, 1500, 1000, 350, 250, 120, 90));
    defineConstant("ep_lengths", c(0, 29800,  1000,  600, 1100,  200,  300,  130,    6,   6,   6,   3,   3));

Category I (`cattle_neutral`) runs the same schedule BACKWARDS in msprime, so it
needs the same numbers in Python. SLiM cannot import a Python module, so this is
a third copy and the three can drift. The guard against that is
``epochs_8_to_12_ticks()``: at Q_scaling = 0.01 it must come out at 2400, which
is the ``handoff_ticks`` every round-3 cattle config states independently. A
transcription slip anywhere in the last five epochs moves it.

Kept free of msprime (and of every other heavy import) so the arithmetic can be
tested without the simulation stack -- helpers/tests stubs those modules out.

WHAT THE Q RESCALING MEANS HERE

`Q` is the config's ``Q_scaling``. SLiM divides every population size and every
epoch duration by it, so at Q = 0.01 the terminal population of 90 real animals
is represented by 9,000 SLiM individuals living 300 ticks instead of 3
generations. Mutation and recombination rates are multiplied by Q to compensate,
which leaves 4*N*mu, 4*N*r and N*s unchanged. Returning the schedule in those
same rescaled units is what lets the coalescent claim to approximate the process
E/F/G actually simulate: a Wright-Fisher population of 9,000 whose entire
membership is sampled, rather than one of 90 sampled 100 times over.

Index 0 of both vectors is dummy padding so that epoch N is at index N, matching
SLiM. Epoch 1 is the burn-in: its size is the ancestral one, and its LENGTH is
unused by the coalescent, which reaches equilibrium on its own rather than
needing a finite number of ticks to get there.
"""

# Epoch:              0      1      2     3     4     5     6     7     8    9    10   11  12
EP_SIZES = (None, 17000, 10000, 7000, 3500, 2500, 2000, 1500, 1000, 350, 250, 120, 90)
EP_LENGTHS = (None, 29800, 1000, 600, 1100, 200, 300, 130, 6, 6, 6, 3, 3)

FIRST_EPOCH = 1
LAST_EPOCH = 12

#: The epoch E/F/G resume from. Everything at or after it is what those arms
#: simulate themselves; everything before it they inherit from the shared
#: `.ep7.ts` checkpoint. Category I simulates the whole thing.
HANDOFF_EPOCH = 8

POPULATION = "Holstein_Friesian"  # the name farm_selection.slim:8 gives it


def _rescaled(value, Q):
    """`value / Q`, the way SLiM writes it (`ep_sizes[i]/Q`)."""
    return value / Q


def epochs_8_to_12_ticks(Q=1.0):
    """Ticks spanned by the epochs E/F/G simulate for themselves.

    At Q = 0.01 this is 600+600+600+300+300 = 2400, the `handoff_ticks` value
    the round-3 cattle configs declare. Asserting it is the only automatic tie
    between this table and the two SLiM copies.
    """
    return sum(_rescaled(EP_LENGTHS[k], Q)
               for k in range(HANDOFF_EPOCH, LAST_EPOCH + 1))


def total_ticks(Q=1.0):
    """Ticks from the start of epoch 2 to the end of epoch 12.

    Epoch 1 is excluded: it is the burn-in, and backwards in time it is the
    ancestral state with no end. 3,354 real generations; 335,400 ticks at
    Q = 0.01.
    """
    return sum(_rescaled(EP_LENGTHS[k], Q)
               for k in range(FIRST_EPOCH + 1, LAST_EPOCH + 1))


def terminal_size(Q=1.0):
    """Population size at sampling: epoch 12, i.e. 9,000 at Q = 0.01."""
    return _rescaled(EP_SIZES[LAST_EPOCH], Q)


def ancestral_size(Q=1.0):
    """Epoch 1's size, in force for all time before epoch 2."""
    return _rescaled(EP_SIZES[FIRST_EPOCH], Q)


def size_changes(Q=1.0):
    """Backwards-time size changes as ``[(time_ago, size, epoch), ...]``.

    Time 0 is the end of epoch 12, which is when SLiM samples. Walking back one
    epoch at a time, the change at ``time_ago`` puts epoch ``k``'s size in force
    for everything older, until the next entry. The terminal size is NOT in the
    list -- it is the population's initial state -- so the first entry is
    epoch 11 and the last is epoch 1, the ancestral size.

    At Q = 0.01 the times come out 300, 600, 1200, 1800, 2400, 15400, 45400,
    65400, 175400, 235400, 335400.
    """
    changes = []
    t = 0.0
    for k in range(LAST_EPOCH - 1, FIRST_EPOCH - 1, -1):
        # Crossing into epoch k means traversing epoch k+1's whole duration.
        t += _rescaled(EP_LENGTHS[k + 1], Q)
        changes.append((t, _rescaled(EP_SIZES[k], Q), k))
    return changes


def describe(Q=1.0):
    """One line per epoch, oldest first, for the stage-1 log."""
    lines = [f'Cattle demography at Q_scaling={Q:g} '
             f'({LAST_EPOCH - FIRST_EPOCH + 1} epochs, {total_ticks(Q):g} ticks '
             f'after the burn-in):']
    lines.append(f'  epoch  {FIRST_EPOCH:>2}  N={ancestral_size(Q):>10,.0f}  '
                 f'(ancestral; extends to the MRCA)')
    for t, size, k in reversed(size_changes(Q)):
        if k == FIRST_EPOCH:
            continue
        lines.append(f'  epoch  {k:>2}  N={size:>10,.0f}  from {t:>10,.0f} ticks ago')
    lines.append(f'  epoch  {LAST_EPOCH:>2}  N={terminal_size(Q):>10,.0f}  '
                 f'from {0:>10,.0f} ticks ago (sampled)')
    return '\n'.join(lines)
