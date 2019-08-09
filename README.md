# HPX_LA
This repository contains a block-based LU factorization implemented in HPX

# Installation
```
git clone https://github.com/jgurhem/HPX_LA.git
cd HPX_LA
mkdir build
cd build
cmake .. -DHPX_DIR=/<path to hpx install dir>/lib64/cmake/HPX
```

# Run
## Use the multithread application
`build/lu --T <nb_blocks> --N <blocksize>`

## Use the distributed application
`mpirun -n <nb_nodes> build/lu_tiled_dist --T <nb_blocks> --N <blocksize> -l <nb_nodes>`


# License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
