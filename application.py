from flask import *
app = Flask(__name__)
import sys
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.indigo_inchi import *
from pymongo import *
import pprint
import base64

INDIGO = Indigo()
RENDERER = IndigoRenderer(INDIGO)
INDIGO_INCHI = IndigoInchi(INDIGO)
DB = Connection('pathway.berkeley.edu', 27334)
CHEMICALS = DB.actv01['chemicals']
REACTIONS = DB.actv01['actfamilies']
FILTER_INFER = json.load(open('../data/infer_ero_pubmed.json'))
FILTER_APPLY = json.load(open('../data/apply_ero_pubmed.json'))
REPORT_REACTIONS = json.load(open('../data/report_reactions.json'))
REPORT_REACTIONS = [(x, y)
                    for x, y in REPORT_REACTIONS if isinstance(y, list) and x != 'Reactions matched']
pr = pprint.PrettyPrinter(indent=2)


def generate_reaction(substrates, products, width=1000, height=300):
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
    INDIGO.setOption('render-image-size', width, height)
    INDIGO.setOption('render-comment', '%s -> %s' %
                     ([x['_id'] for x in substrates], [x['_id'] for x in products]))
    t = RENDERER.renderToBuffer(rxn)
    return 'data:image/png;base64,%s' % base64.b64encode(t)


def generate_ero(query_reaction, width=1000, height=300):
    rxn = INDIGO.loadQueryReaction(query_reaction)
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', width, height)
    t = RENDERER.renderToBuffer(rxn)
    return 'data:image/png;base64,%s' % base64.b64encode(t)


@app.route('/')
def root():
    return render_template('root.html')


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
        filter_apply = FILTER_APPLY.get(rxn_id, [])
        filter_infer = FILTER_INFER.get(rxn_id, [])
    return render_template('rxn.html', reaction=reaction, substrates=substrates, products=products, rxn_img=rxn_img, filter_infer=filter_infer, filter_apply=filter_apply, report_reactions=REPORT_REACTIONS)


if __name__ == '__main__':
    if len(sys.argv) == 2:
        the_file, myport = sys.argv
        app.run(host='0.0.0.0', port=int(myport))
    else:
        app.run(debug=True)  # debug=True will run with reloader enabled
