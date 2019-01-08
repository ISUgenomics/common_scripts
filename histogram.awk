#/usr/bin/awk -f

BEGIN{
    bin_width=100;

}
{
    bin=int(($1-0.0001)/bin_width);
    if( bin in hist){
        hist[bin]+=1
    }else{
        hist[bin]=1
    }
}
END{
    for (h in hist)
        printf " * > %2.2f  ->  %i \n", h*bin_width, hist[h]
}
