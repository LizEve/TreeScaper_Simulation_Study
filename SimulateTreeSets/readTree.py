#!/usr/bin/env python
# Version Oct 10
from __future__ import (division, print_function)
import random
import numpy
# This is currently only compatable with Dendropy 4
import dendropy
import re
from Bio import Phylo
from cStringIO import StringIO

class Node:
	"""
	This class defines the structure of Nodes that will be used to build trees.
	"""
	
	def __init__(self,name="",parent=None,children=None, branchlength = 0):
		"""
		This function initializes the node's:
			- Name
			- Parent
			- Children
			- Branch length
		"""
		self.name = name
		self.brl = branchlength
		if parent is None:
			self.parent = []
		else:
			self.parent = parent
		if children is None:
			self.children = []
		else:
			self.children = children
		
class Tree:
	"""
	Defines a phylogenetic Tree class, consisting of linked Node objects.
	
	Assumes rooted, bifurcating trees.
	"""

	def __init__(self, data, ndict=None, terminal_nodes=None):
		
		# Instantiate a root node
		self.root = Node("root") 
		
		# Create tree from newick string, rooted on self.root
		self.newick_splicer(data.strip(";"), self.root)


	def newick_splicer(self, data, parent):
		"""
		Splices newick tree string (data) to instantiate a Tree object. The string passed
		as data should already have the trailing ";" removed. This tree will be rooted on 
		the node passed as parent.
		"""
		data = data.replace(" ", "")[1: len(data)] 	 #Get rid of all spaces and removes first and last parenthesis
		n = 0
		if data.count(",") != 0: #While there are comma separated taxa
			for key in range(len(data)): #Find the corresponding comma for a given parenthesis (n will be 0 for the correct comma)
				if data[key] == "(":
					n += 1 #Increase index of n by 1 for 1 step into new node
				elif data[key] == ")":
					n -= 1 #Decrease index of n by 1 for 1 step outout node
				elif data[key] == ",":
					if n == 0: #To check for correct comma
						vals = (data[0:key], data[key+1:len(data)-1]) #Break newick into left and right datasets
						for unit in vals: #For each entry of dataset
							if unit.find(":") != -1: #For cases with branch lengths
								d = unit[0:unit.rfind(":")] #get rid of trailing branchlength if provided. rfind from the right side
								node_creater = Node(d, parent = parent) #Create node entry
								node_creater.brl = float(unit[unit.rfind(":")+1:]) #Append branch length of that branch
								parent.children.append(node_creater) #Create children. Hello parent, your children are ...
								self.newick_splicer(d, node_creater) #Recursive function
							else: #For case with no branch lengths
								d=unit
								node_creater = Node(d, parent = parent)
								parent.children.append(node_creater)
								self.newick_splicer(d, node_creater)
						break #Terminate loop, we don't need to look any further


	def print_names(self,node):
		"""
		A method of a Tree object that will print out the names of its terminal nodes. To 
		print out all terminal node names, pass the root as an argument.
		"""
		if node.children == []: #Identifies terminal node
			print (node.name)
		else:
			for child in node.children:
				self.print_names(child)
	
	def list_term_nodes(self,node):
		"""
		A method of a Tree object that will print out the node names and instances for 
		all tips in a tree (or clade) subtended by the provided node. 
		"""
		if node.children == []: # Identifies terminal node
			print(node.name)
			print(node)
		else:	# Internal node
			for child in node.children:
				self.list_term_nodes(child)
	
	
	def inv_edge_len(self,node,edge=0):
		"""
		A method of a Tree object that will return the inverse of the total length from 
		given node to root. 
		"""
		#at root return total
		if node.brl == 0:
			return edge 
		else:
			#add branch length
			edge += 1/float(node.brl)
			return self.inv_edge_len(node.parent,edge)


	def tree_len(self,node):
		"""
		A method to calculate and return total tree length. Pass the root as an argument.
		"""
		tot_len = 0
		if node.children == []: #Terminal branch returns branch length
			return node.brl
		else:
			tot_len += node.brl #Add length of internal branch
			for child in node.children:
				tot_len += self.tree_len(child) #Add length of terminal branch
			return tot_len


	def inv_tree_len(self,node):
		"""
		A method to calculate and return inverse of total tree length.
		"""
		inv_tot_len = 0
		if node.children == []: #Terminal branch returns branch length
			return 1/float(node.brl)
		else:
			if node.brl == 0: #otherwise we get an error for root
				inv_tot_len += 0
				for child in node.children:
					inv_tot_len += self.inv_tree_len(child)
				return inv_tot_len
			else:
				inv_tot_len += 1/float(node.brl) #Add length of internal branch
				for child in node.children:
					inv_tot_len += self.inv_tree_len(child) #Add length of terminal branch
				return inv_tot_len


	def newick(self,node):
		"""
		A method of a Tree object that will print out the Tree as a
		parenthetical string (Newick format).
		"""

		newick = "(" #Opening bracket
		if node.children == []: #Terminal branch returns name
			return node.name + ":" + str(node.brl)
		else:
			for child in node.children:
				if node.children[-1] == child: #Don't add commas to last entry
					newick += self.newick(child)
				else:
					newick += self.newick(child) + "," #Adds commas to non-last entries
			if node.brl != 0:
				newick += "):" + str(node.brl) #Adds closing bracket
			else:
				newick += ")"
			return newick


	def has_grandkids(self,node):
		"""
		Takes a node and will randomy choose a child and return the child node if it has grandchildren
		could add an argument to determine if the child is chosen randomly or based on branch length of children. (for passing shorter branches to NNI moves more often).
		"""
		#if node has children
		if node.children != []:
			#pick random child
			kid = node.children[random.choice([0,1])]
			#if child has children, return grandchildren of node
			if kid.children != []:
				return kid
			else:
				return self.has_grandkids(node)
		#if node doesnt have children, you are at a tip
		else:
			return 0


	def node_dict(self,node,ndict=None):
		"""
		Returns dictionary with all nodes as keys and branch lengths as values.
		"""
		if ndict is None:
			ndict={}
		if node.children == []:		# Terminal node
			ndict[node]=node.brl	# Only terminal branch length returned
			return ndict
		else:						# Internal nodes
			ndict[node]=node.brl 	# Add internal branch length to dictionary
			for child in node.children:
				self.node_dict(child,ndict) # Recursively add all descendant brls
			return ndict

def start_node_exp(tree):
	"""
	Picks node as focus for NNI move based on branch length. Shorter branches will get chosen more often. 
	"""
	# Defines dictionary of nodes and branch lengths for tree
	n_dict = tree.node_dict(tree.root)
	
	# Draws random number from exponential to use when picking node
	goal = numpy.random.exponential(0.1)
	
	# Finds largest brl that's smaller than random exponential (goal)
	# Note for future: could also use smallest brl that's larger than goal
	start_brl = max(brl for brl in list(n_dict.values()) if brl < goal)

	# Finds nodes that matches chosen branch length
	start_nodes_dict = {}
	for key, value in n_dict.items():
		if value == start_brl:
			start_nodes_dict[key]=value
	# Randomly pick node from all nodes that match chosen brl
	start_node = random.choice(start_nodes_dict.keys())		
	return start_node

def start_node_rand(tree):
	"""
	Picks node as focus for NNI move randomly. 
	"""
	# Defines dictionary of nodes and branch lengths for tree
	n_dict = tree.node_dict(tree.root)
	
	# Picks random start node
	start_node = random.choice(n_dict.keys())		
	return start_node

def start_node_filter(tree,node_choice):
	"""
	Serves as wrapper to filter out root and tips.
	"""
	# Start with root
	start_node = tree.root
	
	# Keep picking nodes until the node is not the root or a terminal branch
	while start_node.brl == 0 or start_node.children == []:
		if node_choice == 'exponential':
			#print("exponential")
			start_node = start_node_exp(tree)
		elif node_choice == 'random':
			#print("random")
			start_node = start_node_rand(tree)
		else:
			print("Enter 'exponetial' or 'random' for how you want your starting branch chosen. Check spelling")
	return start_node  

"""
NOTE: Might want to define other functions for picking focal branch for NNI move (e.g. directly proportional to inverse of branch length).
"""

def NNI(orig_tree,node_choice='random'):
	"""
	Does NNI move on random branch, preferentially choosing smaller branches. Returns altered tree. 
	
	Assumes bifurcating tree.
	"""
	tree = Tree(orig_tree.newick(orig_tree.root))
	
	# Store first focal node in c2
	c2 = start_node_filter(tree,node_choice)
	
	# Store 2nd focal node (parent) in p
	p = c2.parent
	
	#print("c2 = "+str(c2.name))
	#print("p = "+str(p.name))
	#assign other child to c1
	
	# Stores other child of p in c1
	for c in p.children:
			if c != c2:
				c1 = c
	
	# c2 will be the node for the brl we choose, so it will always have children. 
	# c1 is the other child of c2's parent
	# Storing branch lengths for c1 and c2. We don't technically need to store these 
	# because they are all still attached to c2, but doing it for clarity for now.
	br1 = c1.brl
	br2 = c2.brl
	
	# Finding children of c2 and storing their branch lengths
	gc1 = c2.children[0]
	gc2 = c2.children[1]
	br3 = gc1.brl
	br4 = gc2.brl
	
	# Remove all children from p
	p.children = []
	
	# Name and instantiate new node
	name = "new_"+str(c2.name)
	new_node = Node(name,parent=p)
	
	# Give it branchlength of c2, then start adding branches
	new_node.brl = br2
	
	# Add new node to parent
	p.children.append(new_node)
	
	# Add c1 to new node
	new_node.children.append(c1)
	c1.parent = new_node
	
	# Reassigning grandkids, one to parent, one to new node. randomly. 
	adopt=random.choice([1,2])
	if adopt == 1:
		p.children.append(gc1)
		gc1.parent = p
		new_node.children.append(gc2)
		gc2.parent = new_node
	elif adopt ==2:
		p.children.append(gc2)
		gc2.parent = p
		new_node.children.append(gc1)
		gc1.parent = new_node

	return tree

def NNI_mult_moves(in_tree,num_moves,node_choice='random',no_dup_start_tree='F', req_min_RF='0', flex_min_RF='F'):
	"""
	Takes an input tree and makes a given number of moves on that tree. outputs a readTree.Tree object.
	Options to not repeat the input tree, and to make a tree a specified normalized RF distance away. If you give a normalized min RF, the number of moves isnt used. 
	I need to redo this function because of the num_moves vs RF issue, but not doing it right now. 
	"""
	new_tree=in_tree
	# If we don't mind getting a duplicate of the start tree
	if float(req_min_RF) == 0:
		if no_dup_start_tree == 'F':
			for move in range(num_moves):
				# Makes given number of moves based on single starting tree.
				new_tree = NNI(new_tree,node_choice)
				dist = rf_unweighted(in_tree,new_tree, normalized='T')[1]

		# If we don't want to ever return to the starting tree. 
		elif no_dup_start_tree == 'T':
			for move in range(num_moves):
				#print("move :"+str(move))
				next_tree = NNI(new_tree,node_choice)
				dist = rf_unweighted(in_tree,next_tree, normalized='T')[1]
				#print("distance :"+str(dist))
				# If new tree is the same as starting tree, make another move on it until it is different. 
				while dist == 0:
					next_tree = NNI(new_tree,node_choice)
					dist = rf_unweighted(in_tree,next_tree, normalized='T')[1]
					#print('redo move')
					#print("distance_x :"+str(dist))
				if dist > 0:
					new_tree = next_tree
					dist = rf_unweighted(in_tree,new_tree, normalized='T')[1]
					#print("new_tree_distance :"+str(dist))
	# If we want the second tree to be a specific RF distance away
	elif float(req_min_RF) != 0:
		req_min_RF=round(req_min_RF,2)
		if flex_min_RF == 'F':
			if no_dup_start_tree == 'F':
				# Do a single NNI move on the new tree
				new_tree = NNI(new_tree,node_choice)
				# Calculate distance between original input tree and new tree
				dist = round(rf_unweighted(in_tree,new_tree, normalized='T')[1],2)
				print("first dist: "+str(dist))
				# While the distance is less than 
				while float(dist) < float(req_min_RF):
					new_tree = NNI(new_tree,node_choice)
					dist = round(rf_unweighted(in_tree,new_tree, normalized='T')[1],2)
					print("distance :"+str(dist))
				while float(dist) > float(req_min_RF):
					new_tree = NNI_mult_moves(in_tree,num_moves,node_choice,no_dup_start_tree, req_min_RF, flex_min_RF)
					dist = round(rf_unweighted(in_tree,new_tree, normalized='T')[1],2)
					print("distance recursion! :"+str(dist))
			elif no_dup_start_tree == 'T':
				print("This option doesn't exist yet. Try no_dup_start_tree == 'T'")
		elif flex_min_RF == 'T':
			if no_dup_start_tree == 'F':
				# Do a single NNI move on the new tree
				new_tree = NNI(new_tree,node_choice)
				# Calculate distance between original input tree and new tree
				dist = rf_unweighted(in_tree,new_tree, normalized='T')[1]
				#print("first dist: "+str(dist))
				# While the distance is less than 
				while float(dist) < float(req_min_RF):
					new_tree = NNI(new_tree,node_choice)
					dist = rf_unweighted(in_tree,new_tree, normalized='T')[1]
					#print("distance :"+str(dist))
			elif no_dup_start_tree == 'T':
				print("This option doesn't exist yet. Try no_dup_start_tree == 'T'")
	return new_tree

def NNI_mult_trees(in_tree,num_out_trees,num_nni_moves,out='file',out_file='outFile.t',node_choice='random',no_dup_start_tree='F'):
	'''
	Does specific number of NNI moves on starting readTree.Tree object and outputs to either nexus file or dendropy tree list object. both with starting tree as first tree in file followed by NNI move trees.
	Does not account for possible duplicates in out file.  
	'''
	# Turns readTree.Tree object into a newick string to read into dendropy
	in_tree_newick = in_tree.newick(in_tree.root)+";"

	# Starts dendropy tree list with starting input tree
	treez = dendropy.TreeList()
	treez.append(dendropy.Tree.get(data=in_tree_newick, schema='newick', rooting='force-rooted'))

	for i in range(num_out_trees):
		# Make NNI move
		new_tree = NNI_mult_moves(in_tree,num_nni_moves,node_choice,no_dup_start_tree)
		# Store tree and read into dendropy
		new_tree = new_tree.newick(new_tree.root)+";"
		nni_tree = dendropy.Tree.get(data=new_tree, schema='newick', rooting='force-rooted')
		# Add to tree list
		treez.append(nni_tree)
	if out == 'file':
		treez.write(path=out_file, schema='nexus')

	elif out == 'list':
		return treez

	else:
		print("Please choose to output of 'list' or 'file'")

'''
In out functions 
'''

def read_nexus(in_file,tree_number=0):
	"""
	Reads in tree from nexus file and returns readTree.Tree object. Picks tree number that you give it. Automatically uses first tree.  
	"""
	# Read in tree from nexus file and transform to newick string
	t = dendropy.Tree.get(path=in_file,schema="nexus",tree_offset=tree_number,rooting='force-rooted')
	# Convert tree to newick string and strip characters.
	u = t.as_string(schema='newick').strip()
	v = u.strip(';')
	w = v.strip('[&R]')
	x = w.strip()
	# Convert to readTree Tree object
	treeOG = Tree(x)
	# Print newick string (not nessecary)
	treeOG.newick(treeOG.root)
	return treeOG

def rand_tree(tips,brl_avg=1,brl_std=None,verbose='T'):
	"""
	Creates random tree to do NNI moves on
	"""
	# Create random tree
	rand_tree = Phylo.BaseTree.Tree.randomized(taxa=tips,branch_length=brl_avg,branch_stdev=brl_std)

	# Convert to newick string, strip trailing whitespace
	rand_newick = rand_tree.format('newick').strip()

	# Remove root branch length of 0
	rand_newick2 = re.sub(':0.00000;',';',rand_newick)

	# Remove node names, use this to read into dendropy if needed
	no_nodes_dp = re.sub('\)n\d+:','):',rand_newick2)

	# Remove trailing ";" use this to read into readTree
	no_nodes = no_nodes_dp.strip(';')
	
	# Print stuff if you want to.
	if verbose == 'T':
		# Print newick string
		print(no_nodes)
		# View tree
		Phylo.draw_ascii(rand_tree)
		
	# Convert to readTree Tree object
	tree_random = Tree(no_nodes)
	return tree_random

def view_phylo(tree_object):
	'''
	Visualize tree with biopython
	'''
	tree_newick = StringIO(tree_object.newick(tree_object.root))
	tree = Phylo.read(tree_newick, "newick")
	print("Tree")
	print(tree_object.newick(tree_object.root))
	Phylo.draw_ascii(tree)

def list_to_out(list_trees,out_file):
	'''
	Takes multiple dendropy tree lists and outputs to a single file.
	'''
	treez = dendropy.TreeList()
	for t in list_trees:
		treez.extend(t)
	treez.write(path=out_file, schema='nexus')

def write_single_tree(in_tree_object,out_file):
	# turn into newick
	newick_tree = in_tree_object.newick(in_tree_object.root) + ";"
	# read into dendropy
	dp_tree = dendropy.Tree.get(data=newick_tree, schema='newick')
	# write to file
	dp_tree.write(path=out_file, schema='nexus')

'''
Tree information functions
'''

def rf_unweighted(tree_object1,tree_object2,normalized='F'): 
	'''
	Gives RF=rf_unweighted[0] and normalized RF=rf_unweighted[1]
	'''
	tree_newick1 = tree_object1.newick(tree_object1.root)+";"
	tree_newick2 = tree_object2.newick(tree_object2.root)+";"
	#print(tree_newick1)
	#print(tree_newick2)
	version = dendropy.__version__.split(".")[0]
	if version == '4':
		taxa = dendropy.TaxonNamespace() #set taxa same for all 
		tree1=dendropy.Tree.get(data=tree_newick1,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
		tree2=dendropy.Tree.get(data=tree_newick2,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
	elif version == '3':
		taxa = dendropy.TaxonSet() #set taxa same for all 
		tree1=dendropy.Tree.get(data=tree_newick1,schema='newick',taxon_set=taxa, rooting='force-rooted')
		tree2=dendropy.Tree.get(data=tree_newick2,schema='newick',taxon_set=taxa, rooting='force-rooted')
	tree1.encode_bipartitions()
	tree2.encode_bipartitions()
	dist=dendropy.calculate.treecompare.symmetric_difference(tree1,tree2)
	if normalized == 'F':
		return dist
	elif normalized == 'T':
		max_RF = 2*(len(taxa)-2)
		norm_dist = dist/max_RF
		both = [dist,norm_dist]
		return both

def rf_weighted(tree_object1,tree_object2): 
	tree_newick1 = tree_object1.newick(tree_object1.root)+";"
	tree_newick2 = tree_object2.newick(tree_object2.root)+";"
	#print(tree_newick1)
	#print(tree_newick2)
	version = dendropy.__version__.split(".")[0]
	if version == '4':
		taxa = dendropy.TaxonNamespace() #set taxa same for all 
		tree1=dendropy.Tree.get(data=tree_newick1,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
		tree2=dendropy.Tree.get(data=tree_newick2,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
	elif version == '3':
		taxa = dendropy.TaxonSet() #set taxa same for all 
		tree1=dendropy.Tree.get(data=tree_newick1,schema='newick',taxon_set=taxa, rooting='force-rooted')
		tree2=dendropy.Tree.get(data=tree_newick2,schema='newick',taxon_set=taxa, rooting='force-rooted')
	tree1.encode_bipartitions()
	tree2.encode_bipartitions()
	dist=dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tree1,tree2)
	return dist

def compare_tree_file(in_file,total_trees,starting_tree_number=0,distance_metric="uRF_norm"):
	'''
	Takes a tree file and compares all trees to first tree in file. 
	For NNI_mult_trees the first tree is the original tree. 
	'''
	if distance_metric == "uRF":
		print("Unweighted Robinson Foulds distance")
	elif distance_metric == "wRF":
		print("Weighted Robinson Foulds distance")
	elif distance_metric == "uRF_norm":
		print("Unweighted Robinson Foulds distance, normalized")
	else:
		print("Input distance metric of choice. 'uRF' for unweighted RF distances, 'wRF' for weighted, 'uRF_norm' for normalized uRF")
	for i in range(0,total_trees):
		in_tree_object=read_nexus(in_file,starting_tree_number)
		new_tree_object=read_nexus(in_file,i)

		if distance_metric == "uRF":
			print("Tree "+str(i)+" : "+str(rf_unweighted(in_tree_object,new_tree_object,'F')))

		elif distance_metric == "uRF_norm":
			un,nm = rf_unweighted(in_tree_object,new_tree_object,'T')
			print("Tree "+str(i)+" : "+str(un)+", "+str(nm))

		elif distance_metric == "RF":
			print("Tree "+str(i)+" : "+str(rf_weighted(in_tree_object,new_tree_object)))
		
def cluster_density_avg(in_file,NNI_trees,starting_tree_number=0):
	'''
	Input number of NNI trees = total trees - start tree
	Compares all NNI created cloud trees to initial seed tree. 
	Counts trees from 0
	'''
	nRF_list = []
	x = starting_tree_number + 1
	y = x + NNI_trees -1
	# print(x,y)
	for i in range(x,y):
		in_tree_object=read_nexus(in_file,starting_tree_number)
		new_tree_object=read_nexus(in_file,i)
		un,nm = rf_unweighted(in_tree_object,new_tree_object,'T')
		nRF_list.append(nm)
	return numpy.mean(nRF_list)


