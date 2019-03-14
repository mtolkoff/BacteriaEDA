import numpy as np
import ast
from flask import request, render_template

from COloNY_functions.FLASK_functions import app

#Produces table of summary information
@app.route('/results_table', methods = ['POST', 'GET'])
def results_table():
    #get summary information
    ratios = request.args.get('ratios')
    ratios = ast.literal_eval(ratios)


    ratios6 = request.args.get('ratios6')
    ratios6 = ast.literal_eval(ratios6)


    ratios7 = request.args.get('ratios7')
    ratios7 = ast.literal_eval(ratios7)


    ratios_other = request.args.get('ratios_other')
    ratios_other = ast.literal_eval(ratios_other)



    x = request.args.get('x_axis')
    x = ast.literal_eval(x)



    y = request.args.get('y_axis')
    y = ast.literal_eval(y)



    OTU = request.args.get('OTU')
    OTU = ast.literal_eval(OTU)



    taxa_name = request.args.get('taxa_name')
    taxa_name = ast.literal_eval(taxa_name.replace("&#39;", "'"))
    taxa_name = ast.literal_eval(taxa_name)
    show_number = min(len(OTU), int(request.args.get('show_number')))
    this_species = OTU.index(int(request.args.get('species')))

    #Measure distance
    species_distance = []

    for i in range(len(OTU)):
        distance_temp = ((x[i] - x[this_species]) ** 2 + (y[i] - y[this_species]) ** 2) ** (.5)
        species_distance.append(distance_temp)

    #sort information to find closest species
    OTU_list = [a for _, a in sorted(zip(species_distance, OTU))]
    ratios_list = [a for _, a in sorted(zip(species_distance, ratios))]
    ratios6_list = [a for _, a in sorted(zip(species_distance, ratios6))]
    ratios7_list = [a for _, a in sorted(zip(species_distance, ratios7))]
    ratios_other_list = [a for _, a in sorted(zip(species_distance, ratios_other))]
    taxa_name_list = [a for _, a in sorted(zip(species_distance, taxa_name))]
    distance_list = np.sort(species_distance)

    
    #Define table. Format so not too many significant digits
    species_table = [['OTU Number', 'Taxa Name', 'Distance', 'Metabolism', 'BRITE Hierarchies', 'Non-BRITE Hierarchies/Unknown Pathway', 'Genetic/Environmental/Cellular Processes']]
#    "{:.8f}".format(float("8.99284722486562e-02"))
    for i in range(show_number):
        temp_row = []
        temp_row.append(OTU_list[i])
        temp_row.append(taxa_name_list[i])
        temp_row.append("{:.1f}".format(float(distance_list[i])))
        temp_row.append(str(ratios_list[i] * 100)[0:4])
        temp_row.append(str(ratios6_list[i] * 100)[0:4])
        temp_row.append(str(ratios7_list[i] * 100)[0:4])
        temp_row.append(str(ratios_other_list[i] * 100)[0:4])

        species_table.append(temp_row)

    return render_template('results_table.html', table = species_table)

