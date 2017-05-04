### VCLab
VCLab (Virtual CALPHAD Laboratory), an open-source software for CALPHAD Calculations written in C++.

###Author
* Xing Wang  <xingwang1991@gmail.com>

###Dependencies

* C++

####Examples

```python
Mode~~~~~~~~~~~~~~~~~~~~~~~= ~~ Equilibrium ~~!\\
Dimension~~~~~~~~~~~~~~~~~= ~~ 1 ~~!\\
Database\_File~~~~~~~~~~~~= ~~ cusnti.TDB ~~! \\
Elements~~~~~~~~~~~~~~~~~~= ~~ cu, sn, ti ~~!\\
Compositions~~~~~~~~~~~~= ~~ 0.1, 0.6, 0.3 ~~!\\
Phases\_Selected~~~~~~~~~~= ~~ all ~~! \\
Phases\_Rejected~~~~~~~~~= ~~ none  ~~! \\
Pressure~~~~~~~~~~~~~~~~~~~= ~~ 101325 ~~!\\
Variables~~~~~~~~~~~~~~~~~~= ~~ Temperature ~~!\\
V\_start~~~~~~~~~~~~~~~~~~~~= ~~ 300 ~~!\\
V\_end~~~~~~~~~~~~~~~~~~~~~= ~~ 2000 ~~!\\
V\_Interval~~~~~~~~~~~~~~~~= ~~ 5 ~~!\\
Global\_Grid\_Interval~~~= ~~ 0.02 ~~!\\
```
Results of phase fractions are stored in phase farctions.txt file. One examples in Cu-Sn-Ti ternary system is:
<img src="documentations/figs/CuSnTi.jpg"/>
