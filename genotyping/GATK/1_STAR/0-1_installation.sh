#0. installation-------------------------------------------------------------------
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git

# Compile
cd STAR/source
make STAR
#make STAR CXXFLAGS_SIMD=sse

# Install
sudo apt-get update
sudo apt-get install rna-star
#----------------------------------------------------------------------------------
