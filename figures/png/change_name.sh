# I have the problem with all of the png file has -1.png in their tail

find . -type f \( -name "*-1.png" -o -name "*-1.PNG" \) -exec sh -c 'mv -i "$1" "${1%-1.png}.png"' _ {} \;
