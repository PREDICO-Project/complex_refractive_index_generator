# complex_refractive_index_generator
Repository to genereate text files with the complex refractive index of different compounds, elements...

The refractive index is complex in the case of X-rays, it is defined as:

$$n = 1 - \delta + i\beta$$

This repository calculates the $\delta$ and $\beta$ values for each energy using the atomic form factors and the mass attenuation coefficient using the following equations:

$$\delta = \frac{\rho N_{A}\lambda^2r_{e}}{M_{a}2\pi}f_{1}$$


$$\beta = \frac{\mu\cdot\lambda}{4\pi}$$


The data inside [form_factors](form_factors) folder have been obtained from [NIST database](https://physics.nist.gov/PhysRefData/FFast/html/form.html). 

* [generate_delta_beta](generate_delta_beta.py) is used to generate the complex refractive index of the elements which atomic form factors are inside the  [form_factors](form_factors) folder.
* [compound_by_chemical_formula](compound_by_chemical_formula.py) enables to calculate the complex refractive index components of a compound using the chemical formula and the density as input.
* [compound_by_weight](compound_by_weight.py) enables to calculate the complex refractive index components of a compound using the weight by mass of the components and the density as input.

The output of all scripts is a text file inside [complex_refractive_index](complex_refractive_index) with the name of the element or compound. The first column of the file is the energy in keV, the second column is the $\delta$ parameter of the refractive index and the third column is the $\beta$ value.

## How to run

Type in the terminal 

```python
python3 generate_delta_beta.py
```
Or:
```python
python3 compound_by_chemical_formula.py
```
Or:
```python
python3 compound_by_weight.py
```
Depending on which script you want to run. **generate_delta_beta.py** and **compound_by_chemical_formula.py** scripts needs to rewrite parameters like the compound name, the density of the compound or the element composition (see the final lines of each script).

