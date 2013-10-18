from flask import *
app = Flask(__name__)
import sys
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.indigo_inchi import *
from pymongo import *

INDIGO = Indigo()
RENDERER = IndigoRenderer(INDIGO)
DB = Connection('pathway.berkeley.edu', 27334)
CHEMICALS = DB.actv01['chemicals']
REACTIONS = DB.actv01['actfamilies']


def generate_reaction(substrates, products):
    output_path = 'static/reaction.png'
    try:
        substrates = [INDIGO.loadMolecule(x) for x in substrates]
        products = [INDIGO.loadMolecule(x) for x in products]
    except IndigoException:
        return False
    rxn = INDIGO.createReaction()
    for s in substrates:
        rxn.addReactant(s)
    for p in products:
        rxn.addProduct(p)
    INDIGO.setOption("render-output-format", "png")
    INDIGO.setOption("render-image-size", 1000, 300)
    RENDERER.renderToFile(rxn, output_path)
    return output_path


def generate_ero(query_reaction):
    output_path = 'static/ero.png'
    rxn = INDIGO.loadQueryReaction(query_reaction)
    INDIGO.setOption("render-output-format", "png")
    INDIGO.setOption("render-image-size", 1000, 300)
    RENDERER.renderToFile(rxn, output_path)
    return output_path


@app.route('/')
def root():
    return render_template('root.html')


@app.route('/rxn/')
@app.route('/rxn/<rxn_id>')
def rxn(rxn_id=None):
    reaction, substrates, products, rxn_img = None, None, None, None
    if rxn_id:
        reaction = REACTIONS.find_one({'_id': long(rxn_id)})
    if reaction:
        products = [CHEMICALS.find_one(product['pubchem'])
                    for product in reaction['enz_summary']['products']]
        product_smiles = [x['SMILES'] for x in products]
        substrates = [CHEMICALS.find_one(product['pubchem'])
                      for product in reaction['enz_summary']['substrates']]
        substrate_smiles = [x['SMILES'] for x in substrates]
        rxn_img = generate_reaction(substrate_smiles, product_smiles)
    return render_template('rxn.html', reaction=reaction, substrates=substrates, products=products, rxn_img=rxn_img)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        the_file, myport = sys.argv
        app.run(host='0.0.0.0', port=int(myport))
    else:
        app.run(debug=True)  # debug=True will run with reloader enabled
