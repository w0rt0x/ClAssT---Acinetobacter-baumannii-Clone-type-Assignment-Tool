# ClAssT - Clone-type Assignment Tool for A.baumannii
ClAssT uses Bloom Filters and Support Vector Machines to classify sequence-reads (.fq-files) or assembled genomes (.fasta or .fna files) to one of the eight international clones of A.baumannii, if the input sequence ist part of any of the eight clones.

The instructions are in the intructions.pdf file.

The requirements are down below.

For security reasons, you need to set up an individual security key and a username/password for a login. More information in the instructions.


# Python Modules - Requirements
## Install Requirements
```
pip install -r requirements.txt
```

## List of used Modules for Python (3.8)
Flask	1.1.2	

Flask-Bcrypt	0.7.1	

Flask-Login	0.5.0	

Flask-WTF	0.14.3	

WTForms	2.3.1	

Werkzeug	1.0.1	

bcrypt	3.1.7	

biopython	1.76	

bitarray	1.2.1	

mmh3	2.5.1	

numpy	1.18.2	

pandas	1.0.3	

requests	2.23.0	

scikit-learn	0.23.1	




# How to run the App: Local Deployment

## MAC/Linux: 

```
$ export FLASK_APP=flaskr

$ export FLASK_ENV=development

$ python Main.py
```

## Windows cmd: 

```
set FLASK_APP=flaskr

set FLASK_ENV=development

python Main.py
```

## More about the Deployment

[Documentation](https://flask.palletsprojects.com/en/master/tutorial/factory/)

# Server Deployment 

[Server Deployment](https://flask.palletsprojects.com/en/master/deploying/)

# Add and remove Filters / Login
```
username: AKEbersberger
password: CavBGmbf20)>
```

