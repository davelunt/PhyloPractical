
from reprophylo import *
import sys

def get_fig_filename(html_name):
    for line in open('.\\phylo_pract_files\\figure.html','r').readlines():
        if line.startswith('<A href'): 
            return line.split('=')[1].split('>')[0]

def cs1_make_clean_alignment():
    # Stating the name of the gene
    pj=Project([Locus('dna','CDS','nd5',['nd5'])])
    
    # Reading the sequence file
    pj.read_denovo(['phylo_pract_files/data/chimp-orang-maca.fasta'],'dna')
    
    # Making sure we are reading the file only once
    if len(pj.records) > 21: 
        pj=Project([Locus('dna','CDS','nd5',['nd5'])])
        pj.read_denovo(['phylo_pract_files/data/chimp-orang-maca.fasta'],'dna')
    print "There are %i sequences in the file"%len(pj.records)
    print
    sys.stdout.flush()
    
    # Say that the sequences are ND5 and get the species and genus
    for r_id in [r.id for r in pj.records]:
        pj.add_feature_to_record(r_id, 'CDS', qualifiers={'gene':'nd5'})
        sample = get_qualifiers_dictionary(pj, "%s_f0"%r_id)['source_original_id'].replace('_',' ')
        genus = sample.split()[0]
        pj.add_qualifier(["%s_f0"%r_id], 'sample',sample)
        pj.add_qualifier(["%s_f0"%r_id], 'genus',genus)
        species = sample
        try:
            int(species[-1])
            species = species[:-1]
        except:
            pass
        species = " ".join(species.split()[:2])
        pj.add_qualifier(["%s_f0"%r_id], 'species',species)
     
       
        
    # Make an alignment
    pj.extract_by_locus()
    print "Running this command to align the sequences:"
    sys.stdout.flush()
    mafft = AlnConf(pj, CDSAlign=False, cmd=pj.defaults['mafft'])
    pj.align([mafft])
    print
    sys.stdout.flush()
    
    # Clean chunks with a lot of missing data
    print "Running this command to remove missing data:"
    sys.stdout.flush()
    trimal = TrimalConf(pj)
    pj.trim([trimal])
    print
    sys.stdout.flush()
    
    # Showing the alignments
    print "\x1b[35mShowing the alignment and the cleaned alignment in new tabs\x1b[0m"
    sys.stdout.flush()
    pj.show_aln('nd5@mafftDefault', id=['sample'])
    pj.show_aln('nd5@mafftDefault@gappyout', id=['sample'])
    return pj
    
def cs1_build_tree(pj):
    raxml = RaxmlConf(pj, threads=1)
    pj.tree([raxml])
    pj.annotate('.\\phylo_pract_files','sample','Pan paniscus2',['sample'],scale = 300, html='.\\phylo_pract_files\\figure.html')
    fig_filename = get_fig_filename('.\\phylo_pract_files\\figure.html')
    if os.path.exists('.\\phylo_pract_files\\primates_tree.png'):
        os.remove('.\\phylo_pract_files\\primates_tree.png')
    os.rename(fig_filename,'.\\phylo_pract_files\\primates_tree.png') 
    return pj, '.\\phylo_pract_files\\primates_tree.png'

def cs1_annotate_tree(pj, species, supports, root_genus):
    # This will correct the root and will reprint the figure
    pj.clear_tree_annotations()
    pj.annotate('.\\phylo_pract_files','genus',root_genus,['sample'],
                leaf_node_color_meta='species', leaf_label_colors=species,
                node_support_dict=supports,
                scale = 300, html='.\\phylo_pract_files\\figure.html')
    fig_filename = get_fig_filename('.\\phylo_pract_files\\figure.html')
    print 'Writing the figure to phylo_pract_files/%s'%'primates_tree.png'
    sys.stdout.flush()
    if os.path.exists('.\\phylo_pract_files\\primates_tree.png'):
        os.remove('.\\phylo_pract_files\\primates_tree.png')
    os.rename(fig_filename,'.\\phylo_pract_files\\primates_tree.png') 
    return '.\\phylo_pract_files\\primates_tree.png'

def cs1_dist_matrix(pj):
    t = pj.ft(pj.trees.keys()[0])
    pant_vs_panp = []
    pongop_vs_pnogoa = []
    leaves = t.get_leaves()
    for i in range(len(leaves)):
        leaf = leaves[i]
        sp = leaf.sample[:-1].replace('pygmaeus pygmaeus','pygmaeus')
        for j in range(i+1,len(leaves)):
            compared_leaf = leaves[j]
            comp_sp = compared_leaf.sample[:-1].replace('pygmaeus pygmaeus','pygmaeus')
            sps = ['Pan paniscus', 'Pan troglodytes']
            if (sp == sps[0] and comp_sp == sps[1]) or (sp == sps[1] and comp_sp == sps[0]):
                pant_vs_panp.append(t.get_distance(leaf.name,compared_leaf.name))
            sps = ['Pongo pygmaeus', 'Pongo abelii']
            if (sp == sps[0] and comp_sp == sps[1]) or (sp == sps[1] and comp_sp == sps[0]):
                pongop_vs_pnogoa.append(t.get_distance(leaf.name,compared_leaf.name))
        
    plt.figure()
    plt.boxplot([pant_vs_panp, pongop_vs_pnogoa])
    plt.xticks([1, 2], ['Within Pan', 'Within Pongo'])
    plt.savefig('.\\phylo_pract_files\\pairwise_tree_distances.png',dpi=300)

def cs2_align_and_tree():
    # Stating the name of the gene
    print 'Reading data...'
    sys.stdout.flush()
    print
    pr=Project([Locus('dna','CDS','coi',['coi'])])
    
    # Reading the sequence file
    pr.read_denovo(['phylo_pract_files/data/AstraptesDNA.fas'],'dna')
    
    # Making sure we are reading the file only once
    if len(pr.records) > 68: 
        pr=Project([Locus('dna','CDS','coi',['coi'])])
        pr.read_denovo(['phylo_pract_files/data/AstraptesDNA.fas'],'dna')
    print "There are %i sequences in the file"%len(pr.records)
    sys.stdout.flush()
    print
    
    # Say that the sequences are coi and get the species, genus, genbank accession and plant
    print 'Sorting metadata...'
    sys.stdout.flush()
    print
    
    plants = {
          'TRIGO' : 'Trigonia' ,
          
          'HIHAMP' : 'Hampea' ,
          
          'YESENN' : 'Senna' ,
          
          'CELT' : 'Celtis' ,
          
          'FABOV' : 'Senna' ,
          
          'INGCUP' : 'Inga or Cupania' ,
          
          'LONCHO' : 'Lonchocarpus' ,
          
          'BYTTNER' : 'Byttneria' ,
          
          'SENNOV' : 'Senna' ,
          
          'LOHAMP' : 'Hampea'
          
          }
    
    for r_id in [r.id for r in pr.records]:
        pr.add_feature_to_record(r_id, 'CDS', qualifiers={'gene':'coi'})
        original_id = get_qualifiers_dictionary(pr, "%s_f0"%r_id)['source_original_id'].replace('\'','')
        accession = original_id.split('_')[0]
        group = original_id.split('_')[-1]
        genus =  original_id.split('_')[1]
        species = ' '.join(original_id.split('_')[1:3])
        pr.add_qualifier(["%s_f0"%r_id], 'original_id',original_id)
        pr.add_qualifier(["%s_f0"%r_id], 'accession',accession)
        pr.add_qualifier(["%s_f0"%r_id], 'group',group)
        pr.add_qualifier(["%s_f0"%r_id], 'plant',plants[group])
        pr.add_qualifier(["%s_f0"%r_id], 'genus',genus)
        pr.add_qualifier(["%s_f0"%r_id], 'species',species)
     
       
        
    # Make an alignment
    print 'Aligning sequences...'
    sys.stdout.flush()
    pr.extract_by_locus()
    print "Running this command to align the sequences:"
    sys.stdout.flush()
    mafft = AlnConf(pr, CDSAlign=False, cmd=pr.defaults['mafft'])
    pr.align([mafft])
    print
    
    # Clean chunks with a lot of missing data
    print 'Cleaning alignment...'
    sys.stdout.flush()
    print "Running this command to remove missing data:"
    sys.stdout.flush()
    trimal = TrimalConf(pr,method_name='gt80', program_name='trimal', cmd='default', trimal_commands={'gt': 0.8})
    pr.trim([trimal])
    print
    
    # Showing the alignments
    #print "\x1b[35mShowing the alignment and the cleaned alignment in new tabs\x1b[0m"
    #pr.show_aln('coi@mafftDefault', id=['original_id'])
    #pr.show_aln('coi@mafftDefault@gt80', id=['original_id'])
    
    print "Reconstructing tree using this command:"
    sys.stdout.flush()
    raxml = RaxmlConf(pr, threads=1)
    pr.tree([raxml])
    pr.clear_tree_annotations()
    pr.annotate('.','mid','mid',['group'],scale = 1000, html='.\\phylo_pract_files\\figure.html')
    fig_filename = get_fig_filename('.\\phylo_pract_files\\figure.html') #Astraptes
    
    if os.path.exists('.\\phylo_pract_files\\astraptes_tree.png'):
        os.remove('.\\phylo_pract_files\\astraptes_tree.png')
    os.rename(fig_filename,'.\\phylo_pract_files\\astraptes_tree.png') 
    return pr, '.\\phylo_pract_files\\astraptes_tree.png'

def cs2_annotate_tree(pr, qual,groups, supports):
    pr.clear_tree_annotations()
    pr.annotate('.\\phylo_pract_files','mid','mid',['group','plant'],
                node_bg_meta=qual, node_bg_color=groups,
                node_support_dict=supports,
                scale = 1000, html='.\\phylo_pract_files\\figure.html')
    fig_filename = get_fig_filename('.\\phylo_pract_files\\figure.html') 
    if os.path.exists('.\\phylo_pract_files\\astraptes_tree.png'):
        os.remove('.\\phylo_pract_files\\astraptes_tree.png')
    os.rename(fig_filename,'.\\phylo_pract_files\\astraptes_tree.png')
    print 'Writing the figure to phylo_pract_files/%s'%'astraptes_tree.png'
    sys.stdout.flush()
    return '.\\phylo_pract_files\\astraptes_tree.png'

def cs3_align_and_tree():
    print "Reading sequence alignmnet..."
    sys.stdout.flush()
    print
    pt = Project([Locus('dna','CDS','ENV',['ENV'])])
    pt.read_alignment('./phylo_pract_files/data/SIVHIV_ENV.fasta','dna','CDS','ENV')
    
    # add strain, host and range of each sample
    print "Sorting metadata..."
    sys.stdout.flush()
    print
    
    pt.add_qualifier_from_source('original_id')
    
    for r in pt.records:
        f = r.features[1]
        orig_id = f.qualifiers['original_id'][0] 
        if 'CP' in orig_id:
            f.qualifiers['strain'] = ['SIVgorGgg']
            f.qualifiers['host'] = ['Gorilla gorilla gorilla']
            f.qualifiers['range'] = ['Cameroon and Gabon']
        elif 'ANT70' == orig_id or 'MVP5180' == orig_id:
            f.qualifiers['strain'] = ['HIV1O']
            f.qualifiers['host'] = ['human']
            f.qualifiers['range'] = ['Cameroon and Gabon']
        elif 'TAN' in orig_id or 'ANT' in orig_id:
            f.qualifiers['strain'] = ['SIVcpzPts']
            f.qualifiers['host'] = ['Pan troglodytes schweinfurthii']
            f.qualifiers['range'] = ['Congo']
        elif 'HXB2' in orig_id or 'U455' in orig_id:
            f.qualifiers['strain'] = ['HIV1M']
            f.qualifiers['host'] = ['human']
            f.qualifiers['range'] = ['Cameroon and Gabon']
        elif 'YBF' in orig_id:
            f.qualifiers['strain'] = ['HIV1N']
            f.qualifiers['host'] = ['human']
            f.qualifiers['range'] = ['Cameroon and Gabon']
        elif 'RBF' in orig_id:
            f.qualifiers['strain'] = ['HIV1P']
            f.qualifiers['host'] = ['human']
            f.qualifiers['range'] = ['Cameroon and Gabon']
        else:
            f.qualifiers['strain'] = ['SIVcpzPtt']
            f.qualifiers['host'] = ['Pan troglodytes troglodytes']
            f.qualifiers['range'] = ['Cameroon and Gabon']
    
    print "Cleaning alignment..."
    sys.stdout.flush()
    trim = TrimalConf(pt)
    pt.trim([trim])
    print
    
    print "Reconstructing tree..."
    sys.stdout.flush()
    raxml = RaxmlConf(pt, threads=1)
    pt.tree([raxml])
    print "Done"
    sys.stdout.flush()
    
    return pt
    
def cs3_annotate_tree(pt, colors):
    t = pt.ft('ENV@ReadDirectly@gappyout@fa')
    ts = TreeStyle()
    ts.mode = "c"
    ts.arc_start = -20
    ts.arc_span = 350
    ts.scale = 200
    ts.show_leaf_name = False
    ts.legend_position=3
    ts.draw_aligned_faces_as_table= False
    
    for c, color in sorted(colors.items()):
        hola = TextFace(" %s "%c)
        hola.background.color = color
        ts.legend.add_face(hola,0)
    
    for s, color in [[' 99-100 ', 'black'],
                     [' 80-99 ', 'gray'],
                     [' 50- 80 ', 'gainsboro']]:
        ts.legend.add_face(TextFace(" "),1)
        ts.legend.add_face(CircleFace(5, color),2)
        ts.legend.add_face(TextFace(" %s "%s),3)
    
    t.ladderize(0)
    t.dist = 0
    t.set_outgroup(t.get_midpoint_outgroup())
    for n in t.traverse():
        ns = NodeStyle()
        ns['size'] = 0
        
        if 50 < n.support < 80:
            n.add_face(CircleFace(4,'gainsboro'),0,position='float')
        if 80 <= n.support < 99:
            n.add_face(CircleFace(4,'gray'),0,position='float')
        if 99 <= n.support < 100:
            n.add_face(CircleFace(4,'black'),0,position='float')
        if n.is_leaf():
            n.add_face(TextFace(n.original_id),0)
            if 'HIV1M' in colors.keys():
                ns['bgcolor'] = colors[n.strain]
            elif 'Congo' in colors.keys():
                ns['bgcolor'] = colors[n.range]
            elif 'human' in colors.keys():
                ns['bgcolor'] = colors[n.host]
        n.set_style(ns)
    
    filename = None
    if 'HIV1M' in colors.keys():
        filename = './phylo_pract_files/HIVSIV_strains.png'
    elif 'Congo' in colors.keys():
        filename = './phylo_pract_files/HIVSIV_range.png'
    elif 'human' in colors.keys():
        filename = './phylo_pract_files/HIVSIV_host.png'
    
    print "Writing the image to %s"%filename
    sys.stdout.flush()
    t.render(filename, tree_style=ts)
    return filename
##