# VEST

VASP Electronic Structure Tool (VEST)   
vest.py is a script to extract data from output file of VASP.

How to use
----
```vim
###########################################################
#         name: VASP Electronic Structure Tool(VEST)      #     
#                                                         #
########### script to extract data from PROCAR ############
# Input file : PROCAR                                     #
#              KPOINTS                                    #
#              POSCAR                                     #
########### script to extract data from EIGENVAL ##########
# VASP version: 5.4.4                                     #
# Input file  : EIGENVAL,(soc,nosoc.megnetic)             #
#               KPOINTS,                                  #
#               POSCAR, fermi.dat                         #
#               KPOINTS.DFT(HSE),                         #
-----------------------------------------------------------
# run command: python3 vest.py                     #
# Author     : Leiwang  updata 2021/04/06                 #
##########################################################
```

 
-----

install
----
1. Firstly, copy vest.py in your bin directory.

```Bash
cp vest.py bin  // Bash   
cd bin    
chmod +x vest.py
```

2. To get the path of your python3

```Bash
which python3               //Bash
/share/pyenv/bin/python3    
```

3. To change the path in the first line of vest.py.

```Vim
#!/share/pyenv/bin/python3  //vim
```

