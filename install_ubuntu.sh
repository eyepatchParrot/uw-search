wget http://releases.llvm.org/5.0.0/clang+llvm-5.0.0-linux-x86_64-ubuntu16.04.tar.xz
tar xvf clang+llvm-5.0.0-linux-x86_64-ubuntu16.04.tar.xz 
mv  clang+llvm-5.0.0-linux-x86_64-ubuntu16.04 ~/clang5

rm clang+llvm-5.0.0-linux-x86_64-ubuntu16.04.tar.xz
sudo apt install -y libomp-dev
sudo apt-get install -y libc++-dev

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x ./Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
rm Miniconda2-latest-Linux-x86_64.sh
conda install jupyter scikit-learn numpy matplotlib pandas

sudo /usr/testbed/bin/mkextrafs /mnt 
