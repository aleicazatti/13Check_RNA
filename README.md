# 13Check_RNA
A tool to re-reference <sup>13</sup>C chemical shifts of RNA

### Installing

1. Download and unzip or clone the repository from https://github.com/BIOS-IMASL/13Check_RNA

2. Execute the following command in your terminal:
 
```
pip install 13Check_RNA/
```

### Using the checkrna module

1. In a python script or notebook shell write:

```
import checkrna
```

2. To check for systematic errors in <sup>13</sup>C chemical shifts of an RNA bmrb entry: 
 
```
checkrna.checkcarbons(bmrb entry number)
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* We thank Jon Wedell, for helping us with the [PyNMRSTAR library](https://github.com/uwbmrb/PyNMRSTAR)
