from flask import *
app = Flask(__name__)
import sys
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.indigo_inchi import *
from pymongo import *
import pprint

INDIGO = Indigo()
RENDERER = IndigoRenderer(INDIGO)
INDIGO_INCHI = IndigoInchi(INDIGO)
DB = Connection('pathway.berkeley.edu', 27334)
CHEMICALS = DB.actv01['chemicals']
REACTIONS = DB.actv01['actfamilies']
FILTER_INFER = json.load(open('../data/infer_ero_pubmed.json'))
FILTER_APPLY = json.load(open('../data/apply_ero_pubmed.json'))
REPORT_REACTIONS = json.load(open('../data/report_reactions.json'))
pr = pprint.PrettyPrinter(indent=2)


def generate_reaction(substrates, products):
    substrate_indigos = [
        INDIGO_INCHI.loadMolecule(x['InChI']) for x in substrates]
    products_indigos = [
        INDIGO_INCHI.loadMolecule(x['InChI']) for x in products]
    rxn = INDIGO.createReaction()
    for s in substrate_indigos:
        rxn.addReactant(s)
    for p in products_indigos:
        rxn.addProduct(p)
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', 1000, 300)
    INDIGO.setOption('render-comment', '%s -> %s' %
                     ([x['_id'] for x in substrates], [x['_id'] for x in products]))
    RENDERER.renderToFile(rxn, 'static/reaction.png')
    return 'reaction.png'


def generate_ero(query_reaction):
    output_path = 'static/ero.png'
    rxn = INDIGO.loadQueryReaction(query_reaction)
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', 1000, 300)
    RENDERER.renderToFile(rxn, output_path)
    return output_path


@app.route('/')
def root():
    return render_template('root.html')


def inchi_to_smiles(inchi):
    return INDIGO_INCHI.loadMolecule(inchi).smiles()


@app.route('/rxn/')
@app.route('/rxn/<rxn_id>')
def rxn(rxn_id=None):
    reaction, substrates, products, rxn_img, filter_apply, filter_infer = None, None, None, None, None, None
    if rxn_id:
        reaction = REACTIONS.find_one({'_id': long(rxn_id)})
    if reaction:
        products = [CHEMICALS.find_one(product['pubchem'])
                    for product in reaction['enz_summary']['products']]
        substrates = [CHEMICALS.find_one(product['pubchem'])
                      for product in reaction['enz_summary']['substrates']]
        rxn_img = generate_reaction(substrates, products)
        filter_apply = pr.pformat(FILTER_APPLY.get(reaction['_id'], []))
        filter_infer = pr.pformat(FILTER_INFER.get(reaction['_id'], []))
    return render_template('rxn.html', reaction=reaction, substrates=substrates, products=products, rxn_img=rxn_img, filter_infer=filter_infer, filter_apply=filter_apply, report_reactions=pr.pformat(REPORT_REACTIONS))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        the_file, myport = sys.argv
        app.run(host='0.0.0.0', port=int(myport))
    else:
        app.run(debug=True)  # debug=True will run with reloader enabled
