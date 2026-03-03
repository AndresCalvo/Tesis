file = "field.dat"

# Data layout: Nt blocks of Nr=101 lines, separated by single blank lines.
# Gnuplot's "index" needs double blank lines, so we use "every" instead.
Nr = 101          # number of r-points per snapshot
Nt = 501          # total number of snapshots (adjust if needed)

print sprintf("Nr = %d, Nt = %d", Nr, Nt)

# --- terminal ---
set term wxt size 1000,700 enhanced
unset output

set xlabel "t"
set ylabel "r"
set zlabel "|A|^2"
set ticslevel 0
set hidden3d
set view 60, 35
unset key

# Uncomment to fix z-range:
# set zrange [0:1]

step = 5
do for [i=0:Nt-1:step] {
    set title sprintf("|A(t,r)|^2  (snapshot %d / %d)", i, Nt-1)
    # "every" skips the blank separator lines automatically for splot;
    # use "index" only with double-blank separated blocks which we don't have.
    # Instead, read the right rows: snapshot i starts at row i*Nr in the data.
    splot file every :::i::i using 1:2:3 with lines
    pause 0.02
}

pause mouse close