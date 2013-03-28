from argparse import ArgumentParser, FileType
from os.path import exists
from sys import argv

import params


def create_parser():
	parser = ArgumentParser()
	parser.add_argument("--params_file",
		type=str,
		metavar="filename",
		default="params.params",
		help="parameters filename")
	parser.add_argument("--job_name",
		type=str,
		help="A general name for the simulation")
	parser.add_argument("--log_file",
		type=str,
		metavar="filename",
		help="log filename")
	parser.add_argument("--pop_size",
		type=int,
		metavar="integer",
		help="population size")
	parser.add_argument( "--U",
		type=float,
		metavar="float",
		help="mutation rate")
	parser.add_argument( "--beta",
		type=float,
		metavar="float",
		help="strain mutation ratio")
	parser.add_argument( "--pi",
		type=int,
		metavar="int",
		help="hypermutaion threshold in fitness")
	parser.add_argument( "--tau",
		type=float,
		metavar="float",
		help="hypermutaion rate fold increase")
	parser.add_argument( "--ticks",
		type=int,
		metavar="integer",
		help="number of ticks")
	parser.add_argument( "--G",
		type=int,
		metavar="integer",
		help="max number of harmful alleles")
	parser.add_argument( "--s",
		type=float,
		metavar="float",
		help="selection coefficient")
	parser.add_argument( "--H",
		type=float,
		metavar="float",
		help="advantage to double mutant")
	parser.add_argument( "--debug",
		action='store_false',
		default=True,
		help="production mode"),
	parser.add_argument( "--console",
		action='store_false',
		default=True,
		help="don't output logging to console")
	parser.add_argument( "--tick_interval",
		type=int,
		metavar="integer",
		help="logging tick interval, 0 for no logging")
	return parser


def parse_args(parser):
	args = parser.parse_args()
	return args


def str2(arg):
	if isinstance(arg, str) or isinstance(arg, unicode):
		return "'"+arg+"'"
	else:
		return str(arg)


def args_and_params():	
	n_args = len(argv)
	if n_args == 2 and not argv[1].startswith('-'):
		parameters = params.load(argv[1])
	else:	
		args = parse_args(create_parser())
		parameters = params.load(args.params_file)
		args = vars(args)
		args = { k: v for k,v in args.items() if v != None }
		parameters.update(args)
	# this is a workaround for sumatra bug (https://groups.google.com/forum/?fromgroups=#!topic/sumatra-users/OIuBWxJF_W0)
	string_to_boolean(parameters,'console')
	string_to_boolean(parameters,'debug')
	return parameters

def string_to_boolean(parameters, field):
	if field in parameters and (isinstance(parameters[field], unicode) or isinstance(parameters[field], str)):
		parameters[field] = parameters[field].lower() == 'true'

if __name__ == '__main__':
	d = args_and_params()
	print params.to_string(d)
