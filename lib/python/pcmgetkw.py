#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vim:syntax=python
#
## @package getkw.py
## @brief getkw -- a simple input parser for Fortran 95
#
# Written by Jonas Juselius <jonas.juselius@chem.uit.no> 
# University of Tromsø, 2006
#
# TODO: 
#       o general cleanup
#       o better python interface
#       o Pretty print inputs back to human readable form
#
# Known bugs: names with '-' mess things up...
#

import sys,os,inspect
import re, string
from copy import deepcopy
from types import *
from pyparsing import \
	Literal, Word, ZeroOrMore, Group, Dict, Optional, removeQuotes, \
	printables, ParseException, restOfLine, alphas, alphanums, nums, \
	pythonStyleComment, oneOf, quotedString, SkipTo, Forward, \
	commaSeparatedList, OneOrMore, Combine, srange, delimitedList, \
	downcaseTokens, line, lineno, StringEnd, Regex

verbose=True
strict=True

class Section:
	"""Section class

	Placehoder for section objects
	"""
	def __init__(self,name,tag=None,req=False, callback=None):
		self.name=name
		self.sect={}
		self.kw={}
		self.tag=tag
		self.req=req
		self.isset=False
		self.callback=callback
		self.fullname=self.name
		if tag != None:
			self.fullname=self.fullname+'<'+self.tag+'>'

	def __cmp__(self, other):
		return cmp(self.name,other.name)

	def __getitem__(self, key):
		if self.sect.has_key(key):
			foo=self.sect
		elif self.kw.has_key(key):
			foo=self.kw
		else:
			return None
		return foo[key]

	def __setitem__(self, k, val):
		if isinstance(val, Section):
			self.sect[k]=val
		elif isinstance(val, Keyword):
			self.kw[k]=val
		else:
			raise TypeError, 'Not a Section or Keyword'
	
	def get(self, k):
		return self.__getitem__(k)

	def set(self, k, val):
		self.__setitem__(k,val)

	def _split_tag(self, key):
		i=string.find(key, '<')
		if i == -1:
			return (key, None)
		j=string.rfind(key, '>')
		if j == -1:
			raise KeyError, 'faulty tag spec'
		return (key[0:i], key[i+1:j])


	def is_set(self, key=None):
		if key is None:
			return self.isset
		if self.kw.has_key(key):
			return self.kw[key].is_set()
		(key, tag)=self._split_tag(key)
		if self.sect.has_key(key):
			sects=self.sect[key]
			if sects.has_key(tag):
				return sects[tag].is_set()
			else:
				raise KeyError, 'missing key'
		return None

	def is_required(self):
		return self.req

	def add_sect(self, sect, set=False):
		s=self.sect
		if s.has_key(sect.name):
			if s[sect.name].has_key(sect.tag):
				print 'Error: Section "%s" already defined!' % \
						(sect.fullname)
				sys.exit(1)
			s[sect.name][sect.tag]=sect
		else:
			s[sect.name]={sect.tag : sect}
		sect.isset=set

	def add_kwkw(self, kw, set=False):
		if not self.kw.has_key(kw.name):
			self.kw[kw.name]=kw
		else:
			print 'Error: Keyword "%s.%s" already defined!' % \
					(self.name, kw.name)
			sys.exit(1)
		kw.isset=set

	def add_kw(self, name, typ, arg=None, req=False, set=False, callback=None):
		if self.kw.has_key(name):
			print 'Error: Keyword "%s.%s" already defined!' % \
					(self.name, kw.name)
			sys.exit(1)
		kw=Keyword(name,typ,arg,req,callback)
		kw.isset=set
		self.kw[name]=kw

	def fetch_kw(self, name):
		if self.kw.has_key(name):
			return self.kw[name]
		return None

	def find_sect(self, path):
		path=string.split(path, '.')
		name=None
		tag=None
		s=self.sect
		pname=''
		for i in path:
			(name, tag)=self._split_tag(i)
			pname=pname+i
			try:
				s=s[name]
			except:
				print 'Invalid section: ', pname
				return None
			else:
				pname=pname+'.'
		return s[tag]

	def get_keyword(self, path):
		path=string.split(path, '.')
		s=self.sect
		k=self.kw
		pname=''
		if len(path) == 1:
			pname=path[0]
			i=pname
		else:
			for i in path[:-1]:
				(name, tag)=self._split_tag(i)
				pname=pname+i
				try:
					s=s[name]
				except:
					str='Invalid section: ' + pname
					raise AttributeError, str
				else:
					pname=pname+'.'
			pname=pname+path[-1:][0]
			i=path[-1:][0]
			k=s[tag].kw
		if k.has_key(i):
			return k[i]
		str='No such key: ' + pname
		raise AttributeError, str

	def getkw(self, path):
		kw=self.get_keyword(path)
		return kw.arg

	def setkw(self, name, arg):
		if self.kw.has_key(name):
			self.kw[name].setkw(arg)
		else:
			print 'Error: invalid kw: ', name

	def fetch_sect(self, name):
		(key, tag)=self._split_tag(name)
		if self.sect.has_key(key):
			if self.sect[key].has_key(tag):
				return self.sect[key][tag]
		return None

	def get_keys(self):
		return self.kw

	def get_sects(self):
		return self.sect

	def status(self):
		return self.isset
	
	def set_status(self, set):
		if set:
			self.isset=True
		else:
			self.isset=False

	def sanitize(self, templ):
		self.equalize(templ)
		self.xvalidate(templ)

	# add missing keys
	def equalize(self, templ):
		for i in templ.kw:
			if not self.kw.has_key(i):
				self.kw[i]=deepcopy(templ.kw[i])
		for i in templ.sect:
			if not self.sect.has_key(i):
				self.sect[i]={None : deepcopy(templ.sect[i][None])}
			for tag in self.sect[i]:
				self.sect[i][tag].equalize(templ.sect[i][None])

	def run_callbacks(self, templ):
		if templ.callback is not None:
			templ.callback(self)
		for i in templ.kw:
			cb=templ.kw[i]
			if cb.callback is not None:
				cb.callback(self.kw[i])
		for i in templ.sect:
			for tag in self.sect[i]:
				self.sect[i][tag].run_callbacks(templ.sect[i][None])

#verify!
	def sanity_check(self,path=None):
		dlm=''
		if path is None:
			path=''
		else:
			path=path+dlm+self.name
			dlm='.'
		if self.req and not self.isset:
			print '>>> Required section not set: %s \n' % (path)
			sys.exit(0)
		for i in self.kw:
			i.sanity_check(path)
		for i in self.sect:
			for j in self.sect[i]:
				j.sanity_check(path)

	#cross-validate against a template
	def xvalidate(self,templ,path=None):
		dlm=''
		if path is None:
			path=''
		else:
			path=path+dlm+self.name
			dlm='.'
		if templ.req and not self.isset:
			print '>>> Required section not set: %s \n' % path
			sys.exit(1)
		for i in self.kw:
			j=templ.fetch_kw(i) 
			if j is None:
				print '>>> Invalid keyword: %s ' % (path+dlm+i)
				sys.exit(1)
			self.kw[i].xvalidate(j,path)
		for i in self.sect:
			j=templ.fetch_sect(i) 
			if j is None:
				print '>>> Invalid section: %s ' % (path+dlm+i)
				sys.exit(1)
			for tag in self.sect[i]:
				self.sect[i][tag].xvalidate(j,path)

	def __str__(self):
		nsect=0
		for i in self.sect:
			for tag in self.sect[i]:
				nsect=nsect+1
		nkw=0
		for i in self.kw:
			nkw=nkw+1
		
		s="SECT %s %d %s\n" % (self.name, nsect, self.isset)
		if self.tag is not None:
			s=s+"TAG T KW %d\n" % (nkw)
			s=s+self.tag+'\n'
		else:
			s=s+"TAG F KW %d\n" % (nkw)

		for i in self.kw:
			s=s+str(self.kw[i])
		for i in self.sect:
			for tag in self.sect[i]:
				s=s+str(self.sect[i][tag])
		return s

	def check_key(self, key):
		(name, tag)=self._split_tag(key)
		try:
			k=self.sect[name][tag]
		except:
			print 'No such key: ', key
			sys.exit(1)
		if k is not None:
			if k.is_set():
				return True
		return False

		
class Keyword:
	"""Placehoder for keyword objects
	"""
	def __init__(self, name, typ, arg=None, req=False, callback=None):
		self.name=name
		self.type=typ
		self.req=req
		self.nargs=None
		self.arg=[]
		self.is_array=False
		self.callback=callback
		
		if re.search('_ARRAY', typ) or typ == 'DATA':
			self.is_array=True
			if arg is None: # unlimited arg length
				self.nargs=-1
			elif isinstance(arg,int): # number of elements in self.arg == arg
				self.nargs=arg  
				arg=None
			if isinstance(arg,tuple) or isinstance(arg,list):
				self.nargs=len(arg)
		else:
			self.nargs=1
		self.setkw(arg)
		self.isset=False # reset the self.isset flag

	def __cmp__(self, other):
		return cmp(self.name,other.name)

	def __getitem__(self, n):
		return self.arg[n]

	def __setitem__(self, n, arg):
		if n > len(self.arg)-1:
			raise IndexError
		self.arg[n]=arg
	
	def setkw(self, arg):
		if arg is None: # init stage
			return

		if self.is_array and self.nargs > 0: 
			if len(arg) != self.nargs:
				print "keyword lenght mismatsh %s(%i): %i" % (self.name,
						self.nargs, len(arg))
				sys.exit(1)
			self.arg=arg
		else:
			self.arg=[arg]
		try:
			self.typecheck()
		except TypeError:
			print 'Invalid argument:', arg
			sys.exit(1)
		self.isset=True

	def get(self):
		if self.is_array:
			return self.arg
		else:
			return self.arg[0]

	def set(self, val):
		self.setkw(val)
	
	def sanity_check(self,path):
		if path is None:
			path=self.name
		else:
			path=path+'.'+self.name
		if self.req  and not self.isset:
			print '>>> Required key not set: %s' % (path)
			if strict:
				sys.exit(1)

	def xvalidate(self,templ,path=None):
		if path is None or path == '':
			path=self.name
		else:
			path=path+'.'+self.name
		if templ.req and not self.isset:
			print '>>> Required key not set: %s' % (path)
			sys.exit(1)
		if templ.type != self.type:
			if self.type == 'INT': 
				self.arg=map(float,self.arg)
				self.type='DBL'
			elif self.type == 'INT_ARRAY':
				self.arg=map(float,self.arg)
				self.type='DBL_ARRAY'
			else:
				print '>>> Invalid data type in: %s' % (path)
				print ' -> wanted %s, got %s' % (templ.type, self.type)
				sys.exit(1)
		if templ.nargs > -1:  #  < 0 == unlimited arg length
			if templ.nargs != len(self.arg):
				print '>>> Invalid data length in: %s' % (path)
				print ' -> wanted %d, got %d' % (templ.nargs, len(self.arg))
				sys.exit(1)
		return True
				
	def typecheck(self):
		if self.arg[0] == 'None':
			return True
		if self.type == None:
			return True
		if (self.type == 'INT' or self.type == 'INT_ARRAY'):
			for i in self.arg:
				if not type(i) == IntType:
					print 'getkw: Not an integer: ', self.name, '=',i
					raise TypeError
		elif (self.type == 'DBL' or self.type == 'DBL_ARRAY'):
			for i in range(len(self.arg)):
				if type(self.arg[i]) == IntType:
					self.arg[i]=float(self.arg[i])
				if not type(self.arg[i]) == FloatType:
					print 'getkw: Not a real: ', self.name, '=', self.arg[i]
					raise TypeError
		elif (self.type == 'BOOL' or self.type == 'BOOL_ARRAY'):
			for i in self.arg:
				if not type(i) == BooleanType:
					print 'getkw: Not a bool: ', self.name, '=',i
					raise TypeError
		elif self.type == 'STR' or self.type == 'STR_ARRAY':
			return True
		elif self.type == 'DATA':
			return True
		else:
			print 'getkw: Unknown type: ', self.name, '=', self.type
			raise TypeError
		return True
	
	def is_type(self, typ):
		if typ == self.type:
			return True
		return False

	def is_set(self):
		return self.isset

	def is_required(self):
		return self.req
	
	def set_status(self, set):
		if set:
			self.isset=True
		else:
			self.isset=False

	def __str__(self):
#        if self.type == 'STR': 
#            print 'foo', self.name, self.arg
		if (self.type == 'STR' or 'STR_ARRAY' or 'DATA') and \
				(self.arg == '' or self.arg == None): # empty string
			nargs=-1 # flags as empty for Fortran code
		else: 
			nargs=len(self.arg)
			tmp=''
			for i in self.arg:
				tmp=tmp+str(i)+'\n'
		s="%s %s %d %s\n" % (self.type, self.name, nargs, self.isset)
		return s+tmp

class Getkw:
	"""Unified interface to sections and keywords.
	Implements a path stack.
	"""
	def __init__(self, top):
		self.top=top
		self.stack=[self.top]
		self.cur=self.stack[0]

	def get_active_section(self):
		return self.cur

	def getkw(self, path):
		return self.cur.getkw(path)

	def get_keyword(self, path):
		retur=self.cur.get_keyword(path)
		return retur

	def setkw(self, path):
		return self.cur.setkw(path)

	def find_sect(self, path):
		return self.cur.find_sect(path)

	def sanitize(self, templ):
		return self.cur.sanitize(templ)

	def run_callbacks(self, templ):
		return self.cur.run_callbacks(templ)

	def push_sect(self, path):
		k=self.cur.find_sect(path)
		if k is None:
			return None
		self.stack.append(k)
		self.cur=self.stack[-1]
		return self.cur

	def pop_sect(self):
		try:
			del self.stack[-1]
		except:
			return None
		self.cur=self.stack[-1]
		return self.cur

	def get_topsect(self):
		return self.top

class GetkwParser:
	"""Implements a class to do the actual parsing of input files and store
	the results in Sections and Keywords. The parseFile() method returns a
	Getkw object.
	"""
	bnf=None
	caseless=False
	yes=re.compile(r'(1|yes|true|on)$',re.I)
	no=re.compile(r'(0|no|false|off)$',re.I)

	def __init__(self,templ=None):
		self.top=Section('toplevel')
		self.stack=[self.top]
		self.cur=self.stack[0]
		self.templ=templ
		self.strg=None
		self.loc=None
		if templ is not None:
			self.path=[self.templ]
		else:
			self.path=None
		if GetkwParser.bnf == None:
			GetkwParser.bnf=self.getkw_bnf()
#        self.parseString=self.bnf.parseString
#        GetkwParser.bnf.setDebug(True)

	def set_caseless(self, arg):
		if arg is True:
			self.caseless=True
		else:
			self.caseless=False

	def parseFile(self,fil):
		self.bnf.parseFile(fil)
		return Getkw(self.top)

	def parseString(self,str):
		self.bnf.parseString(str)
		return Getkw(self.top)

	def add_sect(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		tag=None
		if len(q) > 1:
			if len(q[1]) > 0:
				tag=q[1][0]
		if self.caseless:
			name=name.lower()
		k=Section(name, tag)
		self.cur.add_sect(k, set=True)  
		self.push_sect(k)

	def push_sect(self,k):
		self.stack.append(k)
		self.cur=self.stack[-1]
		if self.templ is not None:
			x=self.path[-1].fetch_sect(k.name)
			if x is None:
				print "Invalid section on line %d: \n%s" % (
						lineno(self.loc,self.strg), line(self.loc,self.strg))
				sys.exit(1)
			self.path.append(x)

	def pop_sect(self,s,l,t):
		if self.templ is not None:  
			del self.path[-1]
		del self.stack[-1]
		self.cur=self.stack[-1]
	
	def store_key(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=q[1]
		if self.caseless:
			name=name.lower()
		if self.templ is None:
			argt=self.fixate_type(arg)
		else:
			k=self.path[-1].fetch_kw(name)
			if k is None:
				print "Unknown keyword '%s' line: %d" % (name, 
						lineno(self.loc,self.strg))
				if strict:
					sys.exit(1)
				argt=None
			else:
				argt=self.check_type(arg,k.type)
		k=Keyword(name,argt,arg)
		self.cur.add_kwkw(k,set=True)

	def store_vector(self,s,l,t):
		q=t.asList()
		self.strg=s
		self.loc=l
		name=q[0]
		arg=q[1:]
		if self.caseless:
			name=name.lower()
		if self.templ is None:
			argt=self.fixate_type(arg[0])
			argt=argt+'_ARRAY'
		else:
			k=self.path[-1].fetch_kw(name)
			if k is None:
				print "Unknown keyword '%s', line: %d" % (name, 
						lineno(self.loc,self.strg))
				if strict:
					sys.exit(1)
				argt=None
			else:
				if k.nargs == -1: 
					pass
				elif len(arg) != k.nargs:
					print "Invalid number of elements for key '%s',\
line: %d" % ( name, lineno(self.loc,self.strg))
					print "  -> %d required, %d given." % (k.nargs, len(arg))
					if strict:
						sys.exit(1)
				argt=self.check_type(arg,k.type)
		k=Keyword(name,argt,arg)
		self.cur.add_kwkw(k,set=True)

	def store_data(self,s,l,t):
		name=t[0]
		self.strg=s
		self.loc=l
		if self.caseless:
			name=name.lower()
		dat=t[1].split('\n')
		dat=dat[:-1] # remove empty element at the end
		arg=[]
		for i in dat:
			arg.append(i.strip())
		k=Keyword(name,'DATA',arg)
		self.cur.add_kwkw(k,set=True)

	def check_type(self, arg, argt):
		print "You hit a bug! Yipeee!"
		sys.exit(1)
		if (isinstance(arg, tuple) or isinstance(arg, list)):
			argt=re.sub('_ARRAY', '', argt)
			for i in arg:
				self.check_type(i, argt)

		if argt == 'INT':
			if not ival.match(arg):
				print 'Invalid type on line %d: Not an int: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt == 'DBL':
			if not dval.match(arg):
				print 'Invalid type on line %d: Not a float: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt == 'BOOL':
			if not lval.match(arg):
				print 'Invalid type on line %d: Not a bool: \n -> %s' % (
						lineno(self.loc,self.strg), line(self.loc,
							self.strg).strip())
				sys.exit(1)
		elif argt != 'STR':
			print 'Invalid type on line %d: Not a %s: \n -> %s' % (
					lineno(self.loc,self.strg), argt, line(self.loc,
						self.strg).strip())
			sys.exit(1)
		return argt

	def fixate_type(self,arg):
		if isinstance(arg, bool):
			return 'BOOL'
		if isinstance(arg, int):
			return 'INT'
		if isinstance(arg, float):
			return 'DBL'
		return 'STR'
	
	def conv_ival(self, s, loc, toks):
		return int(toks[0])

	def conv_dval(self, s, loc, toks):
		return float(toks[0])

	def conv_lval(self, s, loc, toks):
		if GetkwParser.yes.match(toks[0]):
			return True
		elif GetkwParser.no.match(toks[0]):
			return False
		else:
			print "Dynga! Perkele."
		return False

	def getkw_bnf(self):
		sect_begin   = Literal("{").suppress()
		sect_end   = Literal("}").suppress()
		array_begin   = Literal("[").suppress()
		array_end   = Literal("]").suppress()
		tag_begin   = Literal("<").suppress()
		tag_end   = Literal(">").suppress()
		eql   = Literal("=").suppress()
		dmark = Literal('$').suppress()
		end_data=Literal('$end').suppress()
		prtable = alphanums+r'!$%&*+-./<>?@^_|~'
		ival=Regex('[-]?\d+')
		dval=Regex('-?\d+\.\d*([eE]?[+-]?\d+)?')
		lval=Regex('([Yy]es|[Nn]o|[Tt]rue|[Ff]alse|[Oo]n|[Oo]ff)')
	
		# Helper definitions

		kstr= quotedString.setParseAction(removeQuotes) ^ \
				dval ^ ival ^ lval ^ Word(prtable)
		name = Word(alphas+"_",alphanums+"_")
		vec=array_begin+delimitedList(dval ^ ival ^ lval ^ Word(prtable) ^ \
				Literal("\n").suppress() ^ \
				quotedString.setParseAction(removeQuotes))+array_end
		sect=name+sect_begin
		tag_sect=name+Group(tag_begin+name+tag_end)+sect_begin

		# Grammar
		keyword = name + eql + kstr
		vector = name + eql + vec
		data=Combine(dmark+name)+SkipTo(end_data)+end_data
		section=Forward()
		sect_def=(sect | tag_sect ) #| vec_sect)
		input=section | data | vector | keyword 
		section << sect_def+ZeroOrMore(input) + sect_end

		# Parsing actions	
		ival.setParseAction(self.conv_ival)
		dval.setParseAction(self.conv_dval)
		lval.setParseAction(self.conv_lval)
		keyword.setParseAction(self.store_key)
		vector.setParseAction(self.store_vector)
		data.setParseAction(self.store_data)
		sect.setParseAction(self.add_sect)
		tag_sect.setParseAction(self.add_sect)
		sect_end.setParseAction(self.pop_sect)

		bnf=ZeroOrMore(input) + StringEnd().setFailAction(parse_error)
		bnf.ignore(pythonStyleComment)
		return bnf

def parse_error(s,t,d,err):
	print "Parse error, line %d: %s" %  ( lineno(err.loc,err.pstr), 
			line(err.loc,err.pstr))
	sys.exit(1)

######## Convenience routines for callbacks ########

def check_opt(sect,key):
	try:
		k=sect[key]
	except:
		print 'You have a typo in the code for key', key
		sys.exit(1)
	if k is not None:
		if k.is_set():
			return True
	return False

def check_required(list, sect):
	err="Error: Required option '%s' not set in section '%s%s'!" 
	for i in list:
		if not check_opt(sect, i):
			if sect.name == sect.tag or sect.tag is None:
				print err % (i, sect.name, '')
			else:
				print err % (i, sect.name, '<' + sect.tag + '>')
			sys.exit(1)

def check_ignored(list, sect):
	warn="Warning: The '%s' option will be ignored in section '%s%s'." 
	for i in list:
		if check_opt(sect, i):
			if sect.name == sect.tag:
				print warn % (i, sect.name, '')
			else:
				print warn % (i, sect.name, '<' + sect.tag + '>')

####################################################

def test( strng ):
	bnf = GetkwParser()
	try:
		tokens=bnf.parseString(strng)
	except ParseException, err:
		print err.line
		print " "*(err.column-1) + "^"
		print err
	return tokens

if __name__ == '__main__':
	teststr="""
title = foo
string="fooo bar"

defs {
foo=[1,2,3]
bar=99.0
}

defs<apa> { 
foo=[1, 2, 3,
4,5, 6,7,8,9, 
10] 
bar=22.0
}

defs<gorilla> { 
foo=[1, 2, 3,
4,5, 6,7,8,9, 
10] 
bar=22.0
}

verbose=False #(yes|true|on|1)


raboof {
    foo=1
    bar=1

	foobar<gnat>{
		foo=1
		bar=2
	        foobar<gnu>{
		        foo=1
		        bar=2
	        }

	}

    $COORD
      o 0.0 0.0 0.0
      h 1.0 1.0 0.0
      h -1.0 1.0 0.0
    $end
}


"""
#ini = test(teststr)
#print ini.top
#    foo=ini.get_keyword('raboof.foo')
#    print dir(foo)

