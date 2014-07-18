#!/usr/bin/bash

if (( 0 < $# ))
then
    basename=$1
else
    echo "Usage: $0 basename [xy|yz|xz|xyz] [title]"
    exit 1
fi

if (( 1 < $# && 2 <= ${#2} && ${#2} <= 3 ))
then
    space=$2
else
    space="xy"
fi

if (( 2 < $# ))
then
    title=$3
else
    title=$basename
fi

imagedir=$basename$space
mkdir $imagedir 2> /dev/null

for f in frames/frame???.??????.gz
do
    framename=${f:7:12}
    t=$[1${framename:5:3}-1000] # remove leading 0s
    echo $framename
    zcat $f | ./extract_data.rb $space > ${framename}.bin
    sleep 1
    gnuplot -e "imagedir='$imagedir';framename='$framename';title='$title';t=$t" splot${#space}d.gpl
    rm -f ${framename}.bin
done

ffmpeg -y -r 4 -pattern_type glob -i $imagedir'/*.png' -c:v libx264 -vf fps=4 -pix_fmt yuv420p $imagedir.mp4
