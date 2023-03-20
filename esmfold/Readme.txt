### A cookbook for the installation of esmfold ###

Before creating the environment, edit the last line of enviroment_esmfold.yml with the path of your choice. 
To create an environment, use the command : 

conda env create -f enviroment_esmfold.yml

It will create an environment with the required packages and modules. 
Activate it, then :

python -m pip install 'transofmers[esmfold]'

Ignore the warning messages.

And voilà, you can use the program in the created environment.
