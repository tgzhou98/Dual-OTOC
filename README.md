# Dual-OTOC

## Introduction
- `figs`: Figures generated from `Dual-OTOC.nb`
- `julia_src`: Auxiliary `Julia` files support faster calculation on OTOC.
- `Dual-OTOC.nb`: The main file show the numerical result in the paper arXiv:2107.01901. You can open the file with Mathemaitca, and automatically run all initialization cells.

## **Environment Setting**

### In `Dual-OTOC.nb`

#### **Version**
Test on `Mathematica 12.1`

#### **Package `PauliAlgebra`**
This program requires the open source package `PauliAlgebra` from Dr. You-Yi Zhuang(https://github.com/EverettYou/Mathematica-for-physics). For simplicity, this package is put in the root directory. 

#### **Install**
Open Mathemaitca, File -> Install -> install package with path `./PauliAlgebra/PauliAlgebra.m`

#### **Running**
Run the code line by line, or use initial cells evaluations automatically.

### In `julia_src`

#### **Version**
Test on `Julia 1.5.3`

#### **Package**
```
JLD
DelimitedFiles
Plots
LinearAlgebra
LinearAlgebra
Random
Statistics
ProgressMeter
```

## Contributing

Pull requests for new features, bug fixes, and suggestions are welcome!

## License

MIT