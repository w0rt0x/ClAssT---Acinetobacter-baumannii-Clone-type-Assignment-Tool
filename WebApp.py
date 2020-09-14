from flask import Flask, render_template, session, request, redirect, url_for, abort
from Option_field_WebApp import Login
import os
import csv
import json
from search_filter import single_oxa, get_added_genomes, oxa_and_IC_multiprocessing, read_search
from flask_bcrypt import Bcrypt
from flask_login import login_user, login_required, LoginManager
from flask_login import UserMixin
import secrets
import pandas as pd
from Classifier import classify, cut_csv, IC3_classify
import time
from Filter_Editor import add_filter, remove_filter, edit_svm, remove_oxa, add_oxa
import logging
import pickle
# Source Logging and Error Handling
# https://flask.palletsprojects.com/en/1.1.x/logging/
# https://pythonise.com/series/learning-flask/flask-error-handling
# Logging Source: https://stackoverflow.com/questions/17743019/flask-logging-cannot-get-it-to-write-to-a-file
logging.basicConfig(filename='logger.log', level=logging.ERROR)

# init WebApp with flask
app = Flask(__name__)

# reading config file settings
app.config.from_pyfile(r'config/settings.cfg')
bcrypt = Bcrypt(app)
login_manager = LoginManager(app)

# Login System is needed, otherwise the login can be skipped and
# the add/remove filter function is for everyone
#allowed_user = (b'$2b$12$BDm/V.wY54Y/WpllWAJHpumh7J/nnxVgxbm01iomvQPy2v6ixGi2K',
                #b'$2b$12$nB/HaB3etbNhj9FhqayA2eSXuQY4eyh0mbdP/c4aNN7N.WPHmz71i')
with open(r'config/login.txt', 'rb') as fp:
    allowed_user = pickle.load(fp)

# Error Handling:
# https://pythonise.com/series/learning-flask/flask-error-handling


@app.errorhandler(404)
def not_found(e):
    return render_template('404.html')


@app.errorhandler(500)
def not_found(e):
    app.logger.error(f'SERVER ERROR 500 at route {request.url} with error message: {e}')
    app.logger.error(f'Parameters: IC_Lookup{session.get("IC_lookup")}, \n'
                     f'OXA: {session.get("OXA")}, \n'
                     f'QUICK: {session.get("quick")}, \n'
                     f'Filename: {session.get("filename")}, \n'
                     f'Vals OXA: {session.get("vals_oxa")}, \n'
                     f'Vals IC: {session.get("vals_ct")}, \n'
                     f'Hits IC: {session.get("hits_ct")}, \n'
                     f'Time: {session.get("time")}, \n'
                     f'Prediction: {session.get("prediction")}')
    return render_template('500.html')


@app.errorhandler(400)
def not_found(e):
    return render_template('400.html')


@app.errorhandler(401)
def not_found(e):
    return render_template('401.html')


@login_manager.user_loader
def load_user(userid):
    """ validates User """
    # Source: https://gist.github.com/danielfennelly/9a7e9b71c0c38cd124d0862fd93ce217

    if bcrypt.check_password_hash(allowed_user[0], userid):
        user = User()
        user.is_authenticated = True
        user.id = 'User'
        return user
    else:
        return None


class User(UserMixin):
    # Flask-login user class
    # Sources:
    # https://stackoverflow.com/questions/10695093/how-to-implement-user-loader-callback-in-flask-login
    # http://gouthamanbalaraman.com/blog/minimal-flask-login-example.html

    def __init__(self):
        self.id = None
        self._is_authenticated = False

    @property
    def is_authenticated(self):
        return self._is_authenticated

    @is_authenticated.setter
    def is_authenticated(self, val):
        self._is_authenticated = val

    def check_pwd(self, pwd):
        """
        Check user request pwd and update authenticate status.
        """

        if bcrypt.check_password_hash(allowed_user[1], pwd):
            self.is_authenticated = True
        else:
            self.is_authenticated = False


# leads to Tool website
@app.route('/', methods=['GET', 'POST'])
def home():
    """ renders Homepage, gets User parameters and file"""

    # Display added Genomes
    added = get_added_genomes()
    if request.method == 'POST':
        data = request.json

        if data is not None:

            filename = data[-12]
            session['quick'] = data[-11]
            session['IC_lookup'] = data[-10:-1]
            session['OXA'] = data[-1]

            del data[-12:]

            name = r'files/' + str(secrets.token_hex(8)) + filename + '.txt'

            with open(name, 'w') as filehandle:
                for read in data:
                    filehandle.write('%s\n' % read)

            session['filename'] = name

            # Returning a json signal to ajax to redirect to loading page
            # the loading page then triggers the assignment process
            app.logger.info('Assignment started for ' + filename + ', Options: '
                            + str(session.get('IC_lookup', None)) + ', OXA: ' + str(session.get('OXA', None)))
            return json.dumps({'success': True})

        else:
            # Source: https://flask-restplus.readthedocs.io/en/stable/errors.html
            abort(400)

    return render_template('home.html', added=added)


@app.route("/assign")
def assign():
    """ Uses User Options to process the file, returns a signal to the loadingpage to go the the
    result-page when done"""

    # getting user parameters back with session function
    filename = session.get('filename', None)
    IC_lookup = session.get('IC_lookup', None)
    oxa = session.get('OXA', None)
    quick = session.get('quick')

    # user selects cases
    # Case None of the Clonetypes of A.baumannii
    start = time.time()

    if IC_lookup is None or not(os.path.exists(filename)):
        # in case that user types in route of loading screen
        # or file does not exist anymore
        return redirect('/results')

    else:
        # Checking file type
        # if the file is fasta -> concat lines
        ext = filename.split('.')[-2]
        with open(filename) as f:
            reads = f.read().splitlines()

        # Concat Lines if not .fq file
        if ext != 'fq' and ext != 'fastq':
            reads = ''.join(reads)
            reads = reads.split('>')

        else:
            quick = False

        # deleting file
        os.remove(filename)

    if True not in IC_lookup:

        # Oxas also None
        if oxa:

            # Start OXA Assignment, own function
            score_oxa, names_oxa = single_oxa(reads, ext)
            session['vals_oxa'] = score_oxa
            session['names_oxa'] = names_oxa

        else:
            pass
            # Nothing happens

    else:
        # Clonetypes and OXA
        if oxa:
            score_ct, names_ct, hits_ct, score_oxa, names_oxa = oxa_and_IC_multiprocessing(IC_lookup, reads, ext, quick)
            session['vals_oxa'] = score_oxa
            session['names_oxa'] = names_oxa
            session['vals_ct'] = score_ct
            session['names_ct'] = names_ct
            session['hits_ct'] = hits_ct

        else:
            # lookup only for Clonetypes
            # Multiprocessing is only used if the file is a fastq file
            score_ct, names_ct, hits_ct = read_search(IC_lookup, reads, quick)

            # storing values in session for creating plot
            session['vals_ct'] = score_ct
            session['names_ct'] = names_ct
            session['hits_ct'] = hits_ct

        # making prediction
        prediction = classify(r'Training_data/Training_data_IC.csv', score_ct, IC_lookup)
        # Making Label look nicer
        if 'IC' in prediction and len(prediction) == 3:
            prediction = 'International Clonetype ' + prediction[2]

        elif prediction == 'None':
            prediction = 'NONE of the selected Clonetypes or Genomes'

        else:
            pass
        session['prediction'] = prediction

    end = time.time()
    needed = round(end - start, 2)
    session['time'] = str(needed)

    app.logger.info('Assignment done for ' + str(filename) + ', Time needed: ' + str(needed))
    return redirect('/results')

# about page
@app.route('/about')
def about():
    """ returns about page """
    counter = json.load(open(r'filter/OXAs_dict/counter.txt'))
    ids = [*counter]
    r = csv.reader(open(r'Training_data/Training_data_IC.csv'))
    df = pd.DataFrame(data=list(r))
    svm_table = df.to_html(index=False, header=False)
    return render_template('About.html', svm_table=svm_table, oxa_ids=ids)


# add and remove page page
@app.route('/add_and_remove', methods=['GET', 'POST'])
@login_required
def add_and_remove():
    """ returns about page """
    # Pre-OXA-data
    counter = json.load(open(r'filter/OXAs_dict/counter.txt'))
    ids = [*counter]
    if len(ids) > 1:
        allw_oxa_rmv = True
    else:
        allw_oxa_rmv = False

    # SVM Pre-data
    r = csv.reader(open(r'Training_data/Training_data_IC.csv', 'r'))
    svm = list(r)
    header = svm[0]
    svm = svm[1:]

    row_min = len(header[1:-1])
    svm_row = len(svm)
    svm_col = len(svm[0])
    for i in range(len(svm)):
        svm[i] = ','.join(svm[i])
    svm = '\n'.join(svm)

    # Remove filter data
    added = get_added_genomes()

    if added == [None]:
        allow_remove = False
        added = []
    else:
        allow_remove = True

    # Adding lines and Cols for Textarea in 'Add Filter'
    r = csv.reader(open(r'Training_data/Training_data_IC.csv', 'r'))
    svm_add = list(r)
    svm_add = svm_add[1:]

    for i in range(len(svm_add)):
        svm_add[i].insert(-1, 'Score_new')
    names = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7', 'IC8'] + added + ['new']

    for i in range(len(names)):
        names[i] = 'Score_' + names[i]
    new_line = ['Filename'] + names + ['Filtername']

    svm_add.append(new_line)
    svm_add.append(new_line)

    for i in range(len(svm_add)):
        svm_add[i] = ','.join(svm_add[i])
    svm_add = '\n'.join(svm_add)

    # Getting data and executing commands
    if request.method == 'POST':
        data = request.json

        if data[0] == 'SVM':
            # Changing existing filters
            data[1].insert(0, header)
            edit_svm(data[1])
            app.logger.info('SVM data has been edited')

        if data[0] == 'REMOVE':
            # Removing Filter
            remove_filter(data[1])
            app.logger.info('Removed Filter ' + str(data[1]))

        if data[0] == 'ADD':
            # Adding Filter
            header.insert(-1, data[1])
            data[2].insert(0, header)
            add_filter(data[1], data[2], data[3])
            app.logger.info('Added Filter ' + str(data[1]))

        if data[0] == 'REMOVE_OXA':
            # Adding Filter
            app.logger.info('Removing OXA-gene: ' + data[1])
            remove_oxa(data[1])

        if data[0] == 'ADD_OXA':
            app.logger.info('Adding OXA-gene: ' + data[1])
            add_oxa(data[1], data[2])

        # Return JSON Signal to return back to homepage
        return json.dumps({'success': True})

    return render_template('add_and_remove.html',
                           added=added,
                           svm_old=svm,
                           svm_add=svm_add,
                           svm_col=svm_col,
                           svm_row=svm_row,
                           allow_remove=allow_remove,
                           header=header,
                           row_min=row_min,
                           oxa_ids=ids,
                           allow_oxa=allw_oxa_rmv)


# login for add/remove filter
@app.route('/expert_options', methods=['GET', 'POST'])
def expert_options():
    """ returns expert options page """

    error = None
    login_form = Login()
    if login_form.validate_on_submit():
        # Login
        # Source:
        # https://stackoverflow.com/questions/10695093/how-to-implement-user-loader-callback-in-flask-login
        user = User()
        user.id = login_form.name.data
        user.check_pwd(login_form.password.data)

        if user.is_authenticated:
            # Only if Valid Username and pw
            login_user(user)
            app.logger.info('logged in successfully')
            return redirect(url_for('add_and_remove'))

        # error message for invalid Login
        error = 'Invalid Login. This Login is only for a member of the Department for Applied ' \
                'Bioinformatics in Frankfurt'
        app.logger.info('invalid login: User: ' + str(login_form.name.data) + ', PWD: ' + str(login_form.password.data))
        abort(401)

    return render_template('expert_login.html', login_form=login_form, error=error)


@app.route('/results')
def results():
    """ gets Results, creates a Plot and displays them on page with further information"""

    lookup = session.get('IC_lookup')

    # Values of clonetypes, is None if not existing
    values_ct = session.get('vals_ct')
    hits_ct = session.get('hits_ct')
    clonetypes = session.get('names_ct')
    prediction = session.get('prediction')

    # Values of OXAs
    values_oxa = session.get('vals_oxa')
    oxa_names = session.get('names_oxa')

    #Special edge case
    IC3_clade = False
    score_ic3 = 0

    if lookup is not None:

        # validates, if plot for Clonetypes needs to be made or not
        if True not in lookup:
            # if no Clonetypes were selected, no plot will be made
            show_ct = False
            maxi = 0
            prediction = 'N.A.'
            svm_table = None

        else:
            show_ct = True
            maxi = 1

            # Sources for displaying dynamic table:
            # https://stackoverflow.com/questions/52019676/dynamic-table-with-python/52026920
            # https://sarahleejane.github.io/learning/python/2015/08/09/simple-tables-in-webapps-using-flask-and-pandas-with-python.html
            vectors_X, vectors_y = cut_csv(r'Training_data/Training_data_IC.csv', lookup, True)
            df = pd.DataFrame(data=vectors_X)
            svm_table = df.to_html(index=False, header=False)

        # validates, if plot for Oxa-genes needs to be made or not
        if session.get('OXA'):
            show_oxa = True

        else:
            show_oxa = False

        filename = session.get('filename')[22:]
        filename = os.path.splitext(filename)[0]

        return render_template('done.html',
                               show_oxa=show_oxa,  # Display Oxa sector or not
                               results_oxa=values_oxa,
                               oxas=oxa_names,
                               show_CT=show_ct,  # Display Clonetypes sector or not
                               results_ct=values_ct,
                               clonetypes=clonetypes,
                               filename=filename,
                               maxi=maxi,
                               time=session.get('time'),
                               svm_table=svm_table,
                               prediction=prediction)
    else:
        return render_template('done.html',
                               show_oxa=False,
                               results_oxa=values_oxa,
                               oxas='None',
                               maxi_oxa=0,
                               show_CT=False,  # Display Clonetypes sector or not
                               results_ct=0,
                               clonetypes='None',
                               filename='None',
                               maxi=0,
                               time='0')
