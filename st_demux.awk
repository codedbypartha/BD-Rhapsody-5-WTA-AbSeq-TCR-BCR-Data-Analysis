BEGIN {
    for (i=1; i<=length(seq); i++) {
        regexp = regexp sep substr(seq,1,i-1) "." substr(seq,i+1)
        sep = "|"
    }
}
{ rec = rec $0 ORS }
!(NR % 4) {
    if (rec ~ regexp) {
        printf "%s", rec
    }
    rec = ""
}