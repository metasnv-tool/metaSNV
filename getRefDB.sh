#!/bin/bash

download_9() {
    [ -d "db" ] || mkdir -p db
    cd db
    echo "Downloading the archive."
    wget vm-lux.embl.de/~rmuench/files/freeze9.tar.bz2
    echo "This might take up to 10 minutes."
    echo "Extracting files!"
    echo "Patience please..."
    tar xjvf freeze9.tar.bz2
    rm freeze9.tar.bz2
    exit	
}

download_11() {
    [ -d "db" ] || mkdir -p db
    cd db
    echo "Downloading the archive."
    wget vm-lux.embl.de/~rmuench/files/freeze11.tar.bz2
	echo "This might take up to 20 minutes."
    echo "Extracting files!"
    echo "Patience please..."
    tar xjvf freeze11.tar.bz2
    rm freeze11.tar.bz2
    exit
}


echo "You are about to download the metaSNP database"
echo ""
echo "Requirements:"
echo "  - Freeze9:  12G available disc space (1753 genomes)"
echo "  - Freeze11: 33G available disc space (5278 genomes)"
echo ""
echo "Please make a selection by entering the right number:"

select yn in "f9" "f11" "cancel"; do
    case $yn in
        f9 ) download_9 ; break;;
        f11 ) download_11; break;;
        cancel ) exit;;
    esac
done
