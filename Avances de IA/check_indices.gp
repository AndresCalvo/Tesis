file = "field.dat"

# Force a GUI terminal so you still get a window when plotting
set term wxt
unset output

do for [i=0:1000000] {

    # Try to read x (t) from this index block
    stats file index i using 1 nooutput

    if (STATS_records == 0) {
        print sprintf("Index %d: NO numeric x records (empty block or non-numeric like NaN/Inf).", i)
        break
    }

    if (STATS_invalid > 0) {
        print sprintf("Index %d: has %d invalid x values (NaN/Inf etc.).", i, STATS_invalid)
        break
    }

    # Optional: stop after a while if you just want a scan
    if (i > 2000) { print "Scanned 2000 indices, no issue found."; break }
}