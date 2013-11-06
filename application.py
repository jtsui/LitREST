from flask import *
app = Flask(__name__)
import sys
from indigo.indigo import *
from indigo.indigo_renderer import *
from indigo.indigo_inchi import *
from pymongo import *
import pprint
import base64
from collections import defaultdict
import json
import itertools

INDIGO = Indigo()
RENDERER = IndigoRenderer(INDIGO)
INDIGO_INCHI = IndigoInchi(INDIGO)
DB_EROS, EROS, DB, CHEMICALS, REACTIONS, FILTER_INFER, FILTER_APPLY, REPORT_REACTIONS, REACTION_CATEGORIES, REPORT_REACTIONS_SET = None, None, None, None, None, None, None, None, None, None
# these are for development. call initialize() for production data
FILTER_INFER = json.load(open('../data/infer_ero_pubmed.json'))
FILTER_APPLY = json.load(open('../data/apply_ero_pubmed.json'))
REPORT_REACTIONS = json.load(open('../data/report_reactions.json'))
REPORT_REACTIONS = [(x, y)
                    for x, y in REPORT_REACTIONS if isinstance(y, list)]
REACTION_CATEGORIES = defaultdict(list)
REPORT_REACTIONS_SET = {}
all_rxn_ids = set.union(*[set(y) for x, y in REPORT_REACTIONS])
for category, reactions in REPORT_REACTIONS:
    REPORT_REACTIONS_SET[category] = set(reactions)
    REPORT_REACTIONS_SET['^%s' % category] = all_rxn_ids.difference(
        REPORT_REACTIONS_SET[category])
    for reaction in reactions:
        REACTION_CATEGORIES[long(reaction)].append(category)
pr = pprint.PrettyPrinter(indent=2)


def initialize(port, suffix):
    global DB_EROS, EROS, DB, CHEMICALS, REACTIONS, FILTER_INFER, FILTER_APPLY, REPORT_REACTIONS, REACTION_CATEGORIES, REPORT_REACTIONS_SET
    DB_EROS = Connection('pathway.berkeley.edu', 27017)
    EROS = DB_EROS.actv01['eros']
    DB = Connection('pathway.berkeley.edu', port)
    CHEMICALS = DB.actv01['chemicals']
    REACTIONS = DB.actv01['actfamilies']
    FILTER_INFER = json.load(open('../data/infer_ero_pubmed_%s.json' % suffix))
    FILTER_APPLY = json.load(open('../data/apply_ero_pubmed_%s.json' % suffix))
    REPORT_REACTIONS = json.load(
        open('../data/report_reactions_%s.json' % suffix))
    REPORT_REACTIONS = [(x, y)
                        for x, y in REPORT_REACTIONS if isinstance(y, list)]
    REACTION_CATEGORIES = defaultdict(list)
    REPORT_REACTIONS_SET = {}
    all_rxn_ids = set.union(*[set(y) for x, y in REPORT_REACTIONS])
    for category, reactions in REPORT_REACTIONS:
        REPORT_REACTIONS_SET[category] = set(reactions)
        REPORT_REACTIONS_SET['^%s' % category] = all_rxn_ids.difference(
            REPORT_REACTIONS_SET[category])
        for reaction in reactions:
            REACTION_CATEGORIES[long(reaction)].append(category)


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


def generate_ero(query_reaction, width=600, height=200):
    rxn = INDIGO.loadQueryReaction(query_reaction)
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', width, height)
    INDIGO.setOption('render-comment', '')
    t = RENDERER.renderToBuffer(rxn)
    return 'data:image/png;base64,%s' % base64.b64encode(t)


def generate_chem_inchi(inchi, width=150, height=150):
    chem = INDIGO_INCHI.loadMolecule(inchi)
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', width, height)
    INDIGO.setOption('render-comment', '')
    t = RENDERER.renderToBuffer(chem)
    return 'data:image/png;base64,%s' % base64.b64encode(t)


def generate_chems_smiles(smiles, width=450, height=150):
    arr = INDIGO.createArray()
    for x in smiles:
        arr.arrayAdd(INDIGO.loadMolecule(x))
    INDIGO.setOption('render-output-format', 'png')
    INDIGO.setOption('render-image-size', width, height)
    INDIGO.setOption('render-comment', '')
    t = RENDERER.renderGridToBuffer(arr, None, len(smiles))
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
        if filter_apply:
            filter_apply = {'input': filter_apply[0],
                            'inputimg': generate_chem_inchi(filter_apply[0]),
                            'ero': filter_apply[1],
                            'eroimg': generate_ero(EROS.find_one({'_id': filter_apply[1]})['readable'].strip('{').strip('}').strip()),
                            'forward': [(x, generate_chems_smiles(x)) for x in filter_apply[2]['forward'] if x],
                            'reverse': [(x, generate_chems_smiles(x)) for x in filter_apply[2]['reverse'] if x],
                            }
        filter_infer = FILTER_INFER.get(rxn_id, [])
        for res in filter_infer:
            if 'ERO' in res:
                res['EROIMG'] = generate_ero(
                    res['ERO'].strip('{').strip('}').strip())
    return render_template('rxn.html', reaction=reaction, substrates=substrates, products=products, rxn_img=rxn_img, filter_infer=filter_infer, filter_apply=filter_apply, report_reactions=REPORT_REACTIONS, reaction_categories=REACTION_CATEGORIES)


@app.route('/_getrxnids')
def getrxnids():
    checked = json.loads(request.args.get('checked', '', type=str))
    checked_sets = []
    for x in checked:
        checked_sets.append(
            ((x, REPORT_REACTIONS_SET[x]), ('^%s' % x, REPORT_REACTIONS_SET['^%s' % x])))
    result = {}

    if len(checked_sets) == 1:
        result[checked_sets[0][0][0]] = sorted(list(checked_sets[0][0][1]))
        result[checked_sets[0][1][0]] = sorted(list(checked_sets[0][1][1]))
    elif len(checked_sets) > 1:
        for combination in itertools.product(*checked_sets):
            names = [x for x, y in combination]
            intersect = set.intersection(*[y for x, y in combination])
            result[' && '.join(names)] = sorted(list(intersect))
    return render_template('rxnresults.html', result=result)


@app.route('/rxnselect/')
def rxnselect():
    # venn diagram library
    # http://christophermullins.net/page/jquery-venn
    # https://github.com/sidoh/venn
    return render_template('rxnselect.html', categories=sorted(REPORT_REACTIONS_SET.keys()), rxn_ids=sorted(list(set.union(*REPORT_REACTIONS_SET.values()))))


def main():
    if len(sys.argv) == 2:
        the_file, myport, file_suffix, act_port = sys.argv
        initialize(act_port, file_suffix)
        app.run(host='0.0.0.0', port=int(myport))
    elif len(sys.argv) == 1:
        app.run(debug=True)  # debug=True will run with reloader enabled
    else:
        print 'Wrong number of arguments. Usage: python application.py [port] [file suffix] [db port]\nExample: python application.py 27330 Journals10000 27334'
        return

if __name__ == '__main__':
    main()
