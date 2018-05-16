# 13Check_RNA
A tool to evaluate <sup>13</sup>C chemical shifts assignments of RNA

### Installing

&nbsp;&nbsp;&nbsp;&nbsp;0. The Python programming language is a basic requirement. If you don't have it, then install [Anaconda](https://anaconda.org/) or [Miniconda](https://conda.io/miniconda.html). 

1. Download and unzip or clone the repository from https://github.com/BIOS-IMASL/13Check_RNA.

2. Execute the following command in your terminal:
 
```
pip install 13Check_RNA/
```

### Using the checkrna module

[In this Jupyter Notebook](https://github.com/BIOS-IMASL/13Check_RNA/blob/master/examples/ExampleNotebook.ipynb) you can run the following commands and the examples of use, and check-out the output messages and error-warnings. 

0. **In a python script or Jupyter Notebook shell write:**

```python
import checkrna
```

&nbsp;&nbsp;&nbsp;&nbsp;1.0 **To evaluate RNA <sup>13</sup>C chemical shifts from a new nmrstar file not deposited in the [BMRB](http://www.bmrb.wisc.edu/) yet:**

```python
checkrna.checkcarbons('path to the new nmrstar file')
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Where <i>'path to new nmrstar file'</i> is a string.

&nbsp;&nbsp;&nbsp;&nbsp;**Example of use:**

```
checkrna.checkcarbons('example_files/new_nmrstar_file.str')
```

&nbsp;&nbsp;&nbsp;&nbsp;1.1 **To check for systematic errors in <sup>13</sup>C chemical shifts of an RNA BMRB entry:** 

```python
checkrna.checkcarbons(bmrb entry number)
```

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Where <i>bmrb entry number</i> is an integer.

&nbsp;&nbsp;&nbsp;&nbsp;**Example of use:**

```
checkrna.checkcarbons(15869)
```

&nbsp;&nbsp;&nbsp;&nbsp;1.2 **To check for systematic errors in RNA <sup>13</sup>C chemical shifts of a nmrstar file downloaded from the BMRB:** 

```python
checkrna.checkcarbons('path to nmrstar file')
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Where <i>'path to nmrstar file'</i> is a string.


&nbsp;&nbsp;&nbsp;&nbsp;**Example of use:**


```python
checkrna.checkcarbons('example_files/bmr7403.str')
```

### Types of output messages and error-warnings you may expect:

* **If 13C chemical shifts are correct:**


```python
checkrna.checkcarbons(4346)
```

&nbsp;&nbsp;&nbsp;&nbsp;*Nitrogenous base 13C chemical shifts of BMRB id 4346 are correct*\
&nbsp;&nbsp;&nbsp;&nbsp;*Ribose 13C chemical shifts of BMRB id 4346 are correct*

* **If 13C chemical shifts have a systematic error:**

```python
checkrna.checkcarbons(15869)
```
&nbsp;&nbsp;&nbsp;&nbsp;*Nitrogenous base 13C chemical shifts of BMRB id 15869 have a systematic error of -2.69ppm*\
&nbsp;&nbsp;&nbsp;&nbsp;*Ribose 13C chemical shifts of BMRB id 15869 have a systematic error of -2.69ppm*

* **If part of the 13C chemical shifts are correct but another part has a systematic error:**

```python
checkrna.checkcarbons(5932)
```
&nbsp;&nbsp;&nbsp;&nbsp;*Nitrogenous base 13C chemical shifts of BMRB id 5932 are correct*\
&nbsp;&nbsp;&nbsp;&nbsp;*Ribose 13C chemical shifts of BMRB id 5932 have a systematic error of -2.59 ppm*

&nbsp;&nbsp;&nbsp;&nbsp;or

```python
checkrna.checkcarbons(xxxx)
```
&nbsp;&nbsp;&nbsp;&nbsp;*Nitrogenous base 13C chemical shifts of BMRB id xxxx  have a systematic error of -2.71 ppm*\
&nbsp;&nbsp;&nbsp;&nbsp;*Ribose 13C chemical shifts of BMRB id xxxx are correct*

* **If 13C chemical shifts have non-systematic errors (i.e. random errors):**

```python
checkrna.checkcarbons(6239) 
```
&nbsp;&nbsp;&nbsp;&nbsp;*13C chemical shifts of BMRB id 6239 have non-systematic errors*

* **If 13C chemical shifts can not by evaluated because the RNA molecule has not enough reference nuclei:**

```python
checkrna.checkcarbons(xxxx) 
```
&nbsp;&nbsp;&nbsp;&nbsp;*BMRB id xxxx has not enough reference nuclei*

* **If 13C chemical shifts can not by evaluated because the RNA molecule lacks the terminal sequence necessary to apply this method:**

```python
checkrna.checkcarbons(7090)
```
&nbsp;&nbsp;&nbsp;&nbsp;*UserWarning: The method cannot be applied because this RNA molecule has no 5'-GG/3'-C terminal sequence*

* **If BMRB entry doesn't exist:**

```python
checkrna.checkcarbons(919929) 
```
&nbsp;&nbsp;&nbsp;&nbsp;*OSError: Entry '919929' does not exist in the public database.*


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* We thank Jon Wedell, for helping us with the [PyNMRSTAR library](https://github.com/uwbmrb/PyNMRSTAR)



