from flask import Flask, flash, redirect, render_template, request, session, abort
import chem3
from chem3.parser import read_data
from chem3.chemkin import *
import os 

app = Flask(__name__)
 
@app.route("/")
def index():
	return "Flask App!"
 
def get_data(filename):
	a = 1
	b = 2

	test_data_dir = '../../tests/test_data' # os.path.join(os.path.dirname(chem3.__file__), '../tests/test_data')
	db_name = os.path.join(test_data_dir, 'nasa.sqlite')
	file_name = os.path.join(test_data_dir, filename + '.xml')
	
	chemkin_data = chem3.parser.read_data(file_name, db_name)
	system = ReactionSystem(chemkin_data['reactions']['test_mechanism'], chemkin_data['species'])

	data = {}
	data['equations'] = []
	data['species'] = chemkin_data['species']
	for reaction in system.reactions:
		data['equations'].append(reaction.equation)
	return data, system
	
# temporarily use global variables
concs = [2., 1., .5, 1., 1.]
T = 1500
@app.route("/<string:filename>/")
def calculate_rates(filename):
	data, system = get_data(filename)
	rates = system.reaction_rate(concs, T)
	return render_template(
		'test.html', **locals())

 
if __name__ == "__main__":
  print('Start the server')
  app.run(port=5000)

