from carveme import config, project_dir
from carveme.reconstruction.gapfilling import multiGapFill
from carveme.reconstruction.utils import load_media_db
from reframed import load_cbmodel, save_cbmodel
import argparse
import os


def maincall(inputfile, media, mediadb=None, universe=None, universe_file=None, outputfile=None, flavor=None,
         spent=None, verbose=False):

    if verbose:
        print('Loading model...')

    try:
        model = load_cbmodel(inputfile, flavor=flavor)
    except IOError:
        raise IOError('Failed to load model:' + inputfile)

    if spent:
        if verbose:
            print('Loading model for spent medium species...')

        try:
            spent_model = load_cbmodel(spent, flavor=flavor)
        except IOError:
            raise IOError('Failed to load model:' + spent)
    else:
        spent_model = None

    if verbose:
        print('Loading reaction universe...')

    if not universe_file:
        if universe:
            universe_file = "{}{}universe_{}.xml.gz".format(project_dir, config.get('generated', 'folder'), universe)
        else:
            universe_file = project_dir + config.get('generated', 'default_universe')

    try:
        universe_model = load_cbmodel(universe_file, flavor='cobra')
    except IOError:
        if universe:
            raise IOError('Failed to load universe "{0}". Please run build_universe.py --{0}.'.format(universe))
        else:
            raise IOError('Failed to load universe model:' + universe_file)

    if verbose:
        print('Loading media...')

    if not mediadb:
        mediadb = project_dir + config.get('input', 'media_library')

    try:
        media_db = load_media_db(mediadb)
    except IOError:
        raise IOError('Failed to load media database:' + mediadb)

    if verbose:
        m1, n1 = len(model.metabolites), len(model.reactions)
        print('Gap filling for {}...'.format(', '.join(media)))

    max_uptake = config.getint('gapfill', 'max_uptake')
    multiGapFill(model, universe_model, media, media_db, max_uptake=max_uptake, inplace=True,
                 spent_model=spent_model)

    if verbose:
        m2, n2 = len(model.metabolites), len(model.reactions)
        print('Added {} reactions and {} metabolites'.format((n2 - n1), (m2 - m1)))

    if verbose:
        print('Saving SBML file...')

    if not outputfile:
        outputfile = os.path.splitext(inputfile)[0] + '_gapfill.xml'

    if not flavor:
        flavor = config.get('sbml', 'default_flavor')

    save_cbmodel(model, outputfile, flavor=flavor)

    if verbose:
        print('Done.')

def main():
    parser = argparse.ArgumentParser(description="GapFill a metabolic model for a given set of media")

    parser.add_argument('input', metavar='INPUTFILE',
                        help="SBML input file")

    parser.add_argument('-m', '--media', dest='media', required=True, help="List of media (comma-separated)")
    parser.add_argument('--mediadb', help="Media database file")

    parser.add_argument('--spent-medium', metavar='SPECIES', dest='spent',
                        help="Add spent medium compounds generated from given species (SBML model).")

    univ = parser.add_mutually_exclusive_group()
    univ.add_argument('-u', '--universe', dest='universe', help="Pre-built universe model (default: bacteria)")
    univ.add_argument('--universe-file', dest='universe_file', help="Reaction universe file (SBML format)")

    parser.add_argument('-o', '--output', dest='output', type=str, help="SBML output file")

    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', help="Switch to verbose mode")

    sbml = parser.add_mutually_exclusive_group()
    sbml.add_argument('--cobra', action='store_true', help="Input SBML in old cobra format")
    sbml.add_argument('--fbc2', action='store_true', help="Input SBML in sbml-fbc2 format")

    args = parser.parse_args()

    if args.fbc2:
        flavor = 'fbc2'
    elif args.cobra:
        flavor = 'cobra'
    else:
        flavor = config.get('sbml', 'default_flavor')

    maincall(inputfile=args.input,
         media=args.media.split(','),
         mediadb=args.mediadb,
         universe=args.universe,
         universe_file=args.universe_file,
         outputfile=args.output,
         flavor=flavor,
         spent=args.spent,
         verbose=args.verbose)

if __name__ == '__main__':
    main()