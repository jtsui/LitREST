import sys
import base64
import itertools
import parse_utils
from utils import *
from flask import *
from pymongo import *
from indigo.indigo import *
from indigo.indigo_inchi import *
from collections import defaultdict
from indigo.indigo_renderer import *
try:
    import ujson as json
except ImportError:
    import json


app = Flask(__name__)
INDIGO = Indigo()
RENDERER = IndigoRenderer(INDIGO)
INDIGO_INCHI = IndigoInchi(INDIGO)
ERO_DB, DB, CHEM_DB, REACTION_DB = None, None, None, None
FILTER_INFER, FILTER_APPLY = None, None
REPORT_REACTIONS, REACTION_CATEGORIES, REPORT_REACTIONS_SET = None, None, None
HOST = 'pathway.berkeley.edu'
PORT = -1


def load_molecule(smiles_or_inchi):
    if smiles_or_inchi.startswith('InChI='):
        return INDIGO_INCHI.loadMolecule(smiles_or_inchi)
    else:
        return INDIGO.loadMolecule(smiles_or_inchi)


def get_db():
    global DB
    if DB is None:
        DB = Connection(HOST, PORT)
    return DB


def get_ero_db():
    global ERO_DB
    if ERO_DB is None:
        EROS_DB = get_db().actv01['eros']
    return EROS_DB


def get_chem_db():
    global CHEM_DB
    if CHEM_DB is None:
        CHEM_DB = get_db().actv01['chemicals']
    return CHEM_DB


def get_reaction_db():
    global REACTION_DB
    if REACTION_DB is None:
        REACTION_DB = get_db().actv01['actfamilies']
    return REACTION_DB


def initialize(port, suffix):
    global PORT, FILTER_INFER, FILTER_APPLY, REPORT_REACTIONS, REACTION_CATEGORIES, REPORT_REACTIONS_SET
    PORT = port
    FILTER_INFER = json.load(open('../data/filter_infer_%s.json' % suffix))
    FILTER_APPLY = json.load(open('../data/filter_apply_%s.json' % suffix))
    REPORT_REACTIONS = json.load(
        open('../data/filter_report_%s.json' % suffix))
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
    INDIGO.setOption('render-comment', '%s --> %s' %
                     (' + '.join([str(x['_id']) for x in substrates]), ' + '.join([str(x['_id']) for x in products])))
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
    chem = load_molecule(inchi)
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
        reaction = get_reaction_db().find_one({'_id': long(rxn_id)})
        for data in reaction.get('metadata', []):
            sentence = parse_utils.get_sentence(data['sid'])
            for chem in sorted(data['substrates'], key=len):
                sentence = (' %s ' % sentence).replace(
                    ' %s ' % chem, ' <font color="red">%s</font> ' % chem)
            for chem in sorted(data['products'], key=len):
                sentence = (' %s ' % sentence).replace(
                    ' %s ' % chem, ' <font color="blue">%s</font> ' % chem)
            data['sentence'] = sentence.strip()
    if reaction:
        products = [get_chem_db().find_one(product['pubchem'])
                    for product in reaction['enz_summary']['products']]
        for product in products:
            product.update(
                {'img': generate_chem_inchi(product['InChI'], 500, 150)})
        substrates = [get_chem_db().find_one(product['pubchem'])
                      for product in reaction['enz_summary']['substrates']]
        for substrate in substrates:
            substrate.update(
                {'img': generate_chem_inchi(substrate['InChI'], 500, 150)})
        rxn_img = generate_reaction(substrates, products)
        filter_apply = FILTER_APPLY.get(rxn_id, [])
        display_filter_apply = []
        for res in filter_apply:
            display = {}
            display['input'] = res[0]
            display['input_img'] = generate_chems_smiles(res[0])
            display['ero_id'] = res[1]
            ero = get_ero_db().find_one({'_id': res[1]})['readable']
            ero_string = ero.strip('{').strip('}').strip()
            display['ero_img'] = generate_ero(ero_string)
            display['result_img'] = generate_chems_smiles(res[2])
            display['result'] = res[2]
            print res[2]
            display_filter_apply.append(display)
        filter_infer = FILTER_INFER.get(rxn_id, [])
        for res in filter_infer:
            if 'ERO' in res:
                res['EROIMG'] = generate_ero(
                    res['ERO'].strip('{').strip('}').strip())
    return render_template('rxn.html',
                           reaction=reaction,
                           substrates=substrates,
                           products=products,
                           rxn_img=rxn_img,
                           filter_infer=filter_infer,
                           filter_apply=display_filter_apply,
                           report_reactions=REPORT_REACTIONS,
                           reaction_categories=REACTION_CATEGORIES)


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
            names = [x for x, y in combination if not x.startswith('^')]
            names = ' AND '.join(names) or 'Other'
            intersect = set.intersection(*[y for x, y in combination])
            result[names] = sorted(list(intersect))
    return render_template('rxnresults.html', result=result, keys=json.dumps(checked), resultjson=json.dumps(result))


@app.route('/rxnselect/')
def rxnselect():
    return render_template('rxnselect.html', categories=sorted(REPORT_REACTIONS_SET.keys()), rxn_ids=sorted(list(set.union(*REPORT_REACTIONS_SET.values()))))


def main():
    if len(sys.argv) == 4:
        the_file, myport, file_suffix, act_port = sys.argv
        initialize(int(act_port), file_suffix)
        app.run(host='0.0.0.0', port=int(myport), debug=True)
        return
    print 'Wrong number of arguments. Usage: python application.py [port] [file suffix] [db port]'
    print 'Example: python application.py 27330 brenda 27334'

if __name__ == '__main__':
    main()
