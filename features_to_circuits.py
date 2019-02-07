# CRISPR Circuit Generator

# Currently assumes that all TUs (two target sites + gRNA gene) are unique

from sbol import *

SO_GRNA_GENE = 'http://identifiers.org/so/SO:0001264'
SO_NUCLEOTIDE_BINDING_SITE = 'http://identifiers.org/so/SO:0001655'
SO_ENGINEERED_REGION = 'http://identifiers.org/so/SO:0000804'
SO_SGRNA = 'http://identifiers.org/so/SO:0001998'

BIOPAX_RNA_MOLECULE = 'http://www.biopax.org/release/biopax-level3.owl#Rna'

SBO_TEMPLATE = 'http://identifiers.org/biomodels.sbo/SBO:0000645'

SBOL_ORIENTATION = 'http://sbols.org/v2#orientation'

setHomespace('http://hub.sd2e.org/user/sd2e/design')
Config.setOption('sbol_typed_uris', False)
Config.setOption('validate', True)
# Config.setOption('serialization_format', 'rdfxml')

doc = Document()

doc.read('yeast_gates_strains_recursive_fixed_single.xml')

doc2 = Document()

gfp_protein = ComponentDefinition('anno_120009310_protein', BIOPAX_PROTEIN, '1')
gfp_protein.name = 'yeGFP'
gfp_protein.roles = gfp_protein.roles + ['http://purl.obolibrary.org/obo/NCIT_C16586']
doc2.addComponentDefinition(gfp_protein)

gfp_protein.sequence = Sequence('anno_120009310_protein_seq', 'MKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK', SBOL_ENCODING_IUPAC_PROTEIN, '1')

networks = []

for gate in doc.moduleDefinitions:
    print('------------ ' + gate.name)

    network = ModuleDefinition(gate.displayId + '_Network', '1')

    network.name = gate.name + '_Network'
    for role in gate.roles:
        if role != 'http://purl.obolibrary.org/obo/NCIT_C14419':
            network.roles = network.roles + [role]
    networks.append(network)
    
    gate_network = gate.modules.create(network.displayId)
    gate_network.name = network.name
    gate_network.definition = network.identity

    sub_grna_id_to_sub_sites = {}
    sub_grna_id_to_sub_site_orients = {}

    sub_grnas = []

    for gate_plasmid in gate.functionalComponents:
        try:
            plasmid = doc.getComponentDefinition(gate_plasmid.definition)
        except:
            plasmid = None
        
        if plasmid is None:
            print('Plasmid not found: ' + gate_plasmid.definition)
        else:
            seq_annos = plasmid.sequenceAnnotations.getAll()
            seq_annos.sort(key=lambda x: x.locations.getRange().start)

            sub_comps = []
            for seq_anno in seq_annos:
                try:
                    plasmid_comp = plasmid.components.get(seq_anno.component)
                except:
                    plasmid_comp = None

                if plasmid_comp is None:
                    print('Sequence annotation has no component: ' + seq_anno.identity)

                sub_comps.append(plasmid_comp)

            prev_sub_sites = []
            prev_sub_site_orients = []

            prev_revc_sub_grnas = []

            for i in range(0, len(sub_comps)):
                sub_comp = sub_comps[i]

                if sub_comp is not None:
                    try:
                        feat = doc.getComponentDefinition(sub_comp.definition)
                    except:
                        feat = None

                    if feat is not None:
                        print(feat.identity)
                        orient = seq_annos[i].locations.getRange().orientation

                        if SO_GRNA_GENE in feat.roles or feat.identity == 'http://hub.sd2e.org/user/sd2e/design/anno_120009310/1':
                            sub_grnas.append(sub_comp)

                            if orient == SBOL_ORIENTATION_INLINE:
                                sub_grna_id_to_sub_sites[sub_comp.identity] = [] + prev_sub_sites
                                sub_grna_id_to_sub_site_orients[sub_comp.identity] = [] + prev_sub_site_orients
                            elif orient == SBOL_ORIENTATION_REVERSE_COMPLEMENT:
                                sub_grna_id_to_sub_sites[sub_comp.identity] = []
                                sub_grna_id_to_sub_site_orients[sub_comp.identity] = []

                                prev_revc_sub_grnas.append(sub_comp)
                        elif SO_NUCLEOTIDE_BINDING_SITE in feat.roles:
                            prev_sub_sites.append(sub_comp)
                            prev_sub_site_orients.append(orient)

                            for prev_sub_grna in prev_revc_sub_grnas:
                                sub_grna_id_to_sub_sites[prev_sub_grna.identity].append(0, sub_comp)

                            if orient == SBOL_ORIENTATION_INLINE:
                                for prev_sub_grna in prev_revc_sub_grnas:
                                    sub_grna_id_to_sub_site_orients[prev_sub_grna.identity].append(0, SBOL_ORIENTATION_REVERSE_COMPLEMENT)
                            elif orient == SBOL_ORIENTATION_REVERSE_COMPLEMENT:
                                for prev_sub_grna in prev_revc_sub_grnas:
                                    sub_grna_id_to_sub_site_orients[prev_sub_grna.identity].append(0, SBOL_ORIENTATION_INLINE)

                        elif SO_TERMINATOR in feat.roles:
                            if orient == SBOL_ORIENTATION_INLINE:
                                prev_sub_sites.clear()
                                prev_sub_site_orients.clear()
                            elif orient == SBOL_ORIENTATION_REVERSE_COMPLEMENT:
                                prev_revc_sub_grnas.clear()

    seq_to_network_grna_molecule = {}
    seq_to_network_grna_tus = {}

    network_grna_tus = []

    for sub_grna in sub_grnas:
        grna = doc.getComponentDefinition(sub_grna.definition)

        # Create TU ComponentDefinition and FunctionalComponent

        grna_tu_id = grna.displayId + '_tu'

        n = 0
        id_is_unique = False
        while not id_is_unique:
            try:
                if n > 0:
                    grna_tu = ComponentDefinition('_'.join([grna_tu_id, repr(n)]), BIOPAX_DNA, '1')
                else:
                    grna_tu = ComponentDefinition(grna_tu_id, BIOPAX_DNA, '1')
                doc2.addComponentDefinition(grna_tu)

                id_is_unique = True
            except:
                n = n + 1
                
        grna_tu.name = grna.name.replace('Gene', 'TU')
        grna_tu.roles = grna_tu.roles + [SO_ENGINEERED_REGION]

        network_grna_tu = network.functionalComponents.create(grna_tu.displayId)

        network_grna_tu.name = grna_tu.name
        network_grna_tu.definition = grna_tu.identity

        network_grna_tus.append(network_grna_tu)

        sub_sites = sub_grna_id_to_sub_sites[sub_grna.identity]
        sub_site_orients = sub_grna_id_to_sub_site_orients[sub_grna.identity]

        prev_tu_site = None

        for i in range(0, len(sub_sites)):
            site = doc.getComponentDefinition(sub_sites[i].definition)

            site_seq = site.sequence.elements.upper()

            if site_seq not in seq_to_network_grna_tus:
                seq_to_network_grna_tus[site_seq] = []
            seq_to_network_grna_tus[site_seq].append(network_grna_tu)

            n = 0
            id_is_unique = False
            while not id_is_unique:
                try:
                    if n > 0:
                        tu_site = grna_tu.components.create('_'.join([site.displayId, repr(n)]))
                    else:
                        tu_site = grna_tu.components.create(site.displayId)

                    id_is_unique = True
                except:
                    n = n + 1
            
            tu_site.name = site.name
            tu_site.definition = site.identity

            tu_site_annotation = grna_tu.sequenceAnnotations.create(tu_site.displayId + '_annotation')

            tu_site_annotation.name = tu_site.name + ' Annotation'
            tu_site_annotation.component = tu_site.identity
            
            tu_site_location = tu_site_annotation.locations.createGenericLocation(tu_site.displayId + '_location')

            tu_site_location.setPropertyValue(SBOL_ORIENTATION, sub_site_orients[i])

            if prev_tu_site is not None:
                tu_site_constraint = grna_tu.sequenceConstraints.create('_precedes_'.join([prev_tu_site.displayId, tu_site.displayId]))

                tu_site_constraint.name = ' Precedes '.join([prev_tu_site.name, tu_site.name])
                tu_site_constraint.subject = prev_tu_site.identity
                tu_site_constraint.object = tu_site.identity
                tu_site_constraint.restriction = SBOL_RESTRICTION_PRECEDES

            prev_tu_site = tu_site

        n = 0
        id_is_unique = False
        while not id_is_unique:
            try:
                if n > 0:
                    tu_gene = grna_tu.components.create('_'.join([grna.displayId, repr(n)]))
                else:
                    tu_gene = grna_tu.components.create(grna.displayId)

                id_is_unique = True
            except:
                n = n + 1

        tu_gene.name = grna.name
        tu_gene.definition = grna.identity

        tu_gene_annotation = grna_tu.sequenceAnnotations.create(tu_gene.displayId + '_annotation')

        tu_gene_annotation.name = tu_gene.name + ' Annotation'
        tu_gene_annotation.component = tu_gene.identity
            
        tu_gene_annotation.locations.createGenericLocation(tu_gene.displayId + '_location')

        if prev_tu_site is not None:
            tu_gene_constraint = grna_tu.sequenceConstraints.create('_precedes_'.join([prev_tu_site.displayId, tu_gene.displayId]))

            tu_gene_constraint.name = ' Precedes '.join([prev_tu_site.name, tu_gene.name])
            tu_gene_constraint.subject = prev_tu_site.identity
            tu_gene_constraint.object = tu_gene.identity
            tu_gene_constraint.restriction = SBOL_RESTRICTION_PRECEDES

        # Create gRNA/GFP ComponentDefinition and FunctionalComponent if needed

        if grna.identity == 'http://hub.sd2e.org/user/sd2e/design/anno_120009310/1':
            network_grna_molecule = network.functionalComponents.create(gfp_protein.displayId)

            network_grna_molecule.name = gfp_protein.name
            network_grna_molecule.definition = gfp_protein.identity

            grna_molecule = gfp_protein
        else:
            grna_seq = grna.sequence.elements.upper()

            grna_molecule_id = grna.displayId + '_rna'

            try:
                network_grna_molecule = network.functionalComponents.get(grna_molecule_id)
                # grna_molecule = doc2.getComponentDefinition('/'.join([getHomespace(), grna_molecule_id, '1']))
            except:
                grna_molecule = ComponentDefinition(grna_molecule_id, BIOPAX_RNA_MOLECULE, '1')

                doc2.addComponentDefinition(grna_molecule)

                grna_molecule.name = grna.name.replace(' Gene', '')
                grna_molecule.roles = grna_molecule.roles + [SO_SGRNA]

                network_grna_molecule = network.functionalComponents.create(grna_molecule.displayId)

                network_grna_molecule.name = grna_molecule.name
                network_grna_molecule.definition = grna_molecule.identity

                seq_to_network_grna_molecule[grna_seq] = network_grna_molecule

                grna_molecule.sequence = Sequence(grna_molecule.displayId + '_seq', grna_seq.replace('T', 'U'), SBOL_ENCODING_IUPAC, '1')
                grna_molecule.sequence.name = grna_molecule.name + ' Sequence'

        # Create network interaction for production of gRNA/GFP molecule from TU

        grna_tu_production = network.interactions.create('_produces_'.join([network_grna_tu.displayId, network_grna_molecule.displayId]))
        grna_tu_production.name = ' Produces '.join([grna_tu.name, grna_molecule.name])
        grna_tu_production.types = [SBO_GENETIC_PRODUCTION]

        grna_tu_template = grna_tu_production.participations.create(network_grna_tu.displayId)
        grna_tu_template.name = network_grna_tu.name
        grna_tu_template.participant = network_grna_tu.identity
        grna_tu_template.roles = grna_tu_template.roles + [SBO_TEMPLATE]

        grna_molecule_product = grna_tu_production.participations.create(network_grna_molecule.displayId)
        grna_molecule_product.name = network_grna_molecule.name
        grna_molecule_product.participant = network_grna_molecule.identity
        grna_molecule_product.roles = grna_molecule_product.roles + [SBO_PRODUCT]

    # Create network interactions for repression of TUs by gRNA molecules with corresponding DNA sequences

    for seq in seq_to_network_grna_molecule:
        network_grna_molecule = seq_to_network_grna_molecule[seq]

        for network_grna_tu in seq_to_network_grna_tus[seq]:
            grna_molecule = doc2.getComponentDefinition(network_grna_molecule.definition)
            grna_tu = doc2.getComponentDefinition(network_grna_tu.definition)

            grna_tu_repression = network.interactions.create('_represses_'.join([network_grna_molecule.displayId, network_grna_tu.displayId]))
            grna_tu_repression.name = ' Represses '.join([grna_molecule.name, grna_tu.name])
            grna_tu_repression.types = [SBO_INHIBITION]

            grna_molecule_repressor = grna_tu_repression.participations.create(network_grna_molecule.displayId)
            grna_molecule_repressor.name = network_grna_molecule.name
            grna_molecule_repressor.roles = grna_molecule_repressor.roles + [SBO_INHIBITOR]
            grna_molecule_repressor.participant = network_grna_molecule.identity

            grna_tu_repressed = grna_tu_repression.participations.create(network_grna_tu.displayId)
            grna_tu_repressed.name = network_grna_tu.name
            grna_tu_repressed.roles = grna_tu_repressed.roles + [SBO_INHIBITED]
            grna_tu_repressed.participant = network_grna_tu.identity

    for i in range(0, len(network_grna_tus)):
        network_grna_tu = network_grna_tus[i]
        grna_tu = doc2.getComponentDefinition(network_grna_tu.definition)

        sub_grna = sub_grnas[i]
        grna = doc.getComponentDefinition(sub_grna.definition)

        # Create TU interactions ModuleDefinition and Module

        grna_tu_interactions = ModuleDefinition(network_grna_tu.displayId + '_interactions', '1')
        doc2.addModuleDefinition(grna_tu_interactions)

        grna_tu_interactions.name = grna_tu.name + ' Interactions'

        network_grna_tu_interactions = network.modules.create(grna_tu_interactions.displayId)

        network_grna_tu_interactions.name = grna_tu_interactions.name
        network_grna_tu_interactions.definition = grna_tu_interactions.identity

        interactions_grna = grna_tu_interactions.functionalComponents.create(grna.displayId)

        interactions_grna.name = grna.name
        interactions_grna.definition = grna.identity

        if grna.identity == 'http://hub.sd2e.org/user/sd2e/design/anno_120009310/1':
            interactions_out_molecule = grna_tu_interactions.functionalComponents.create(gfp_protein.displayId)

            interactions_out_molecule.name = gfp_protein.name
            interactions_out_molecule.definition = gfp_protein.identity

            out_molecule = gfp_protein
        else:
            grna_seq = grna.sequence.elements.upper()

            network_grna_molecule = seq_to_network_grna_molecule[grna_seq]
            out_molecule = doc2.getComponentDefinition(network_grna_molecule.definition)

            interactions_out_molecule = grna_tu_interactions.functionalComponents.create(out_molecule.displayId)

            interactions_out_molecule.name = out_molecule.name
            interactions_out_molecule.definition = out_molecule.identity

        interactions_production = grna_tu_interactions.interactions.create('_produces_'.join([interactions_grna.displayId, interactions_out_molecule.displayId]))
        interactions_production.name = ' Produces '.join([grna.name, out_molecule.name])
        interactions_production.types = [SBO_GENETIC_PRODUCTION]

        interactions_template = interactions_production.participations.create(interactions_grna.displayId)
        interactions_template.name = interactions_grna.name
        interactions_template.participant = interactions_grna.identity
        interactions_template.roles = interactions_template.roles + [SBO_TEMPLATE]

        interactions_product = interactions_production.participations.create(interactions_out_molecule.displayId)
        interactions_product.name = interactions_out_molecule.name
        interactions_product.participant = interactions_out_molecule.identity
        interactions_product.roles = interactions_product.roles + [SBO_PRODUCT]
         
        sub_sites = sub_grna_id_to_sub_sites[sub_grna.identity]

        for sub_site in sub_sites:
            site = doc.getComponentDefinition(sub_site.definition)
            site_seq = site.sequence.elements.upper()

            if site_seq in seq_to_network_grna_molecule:
                network_grna_molecule = seq_to_network_grna_molecule[site_seq]
                in_molecule = doc2.getComponentDefinition(network_grna_molecule.definition)

                try:
                    interactions_in_molecule = grna_tu_interactions.functionalComponents.get(in_molecule.displayId)
                except:
                    interactions_in_molecule = grna_tu_interactions.functionalComponents.create(in_molecule.displayId)

                interactions_in_molecule.name = in_molecule.name
                interactions_in_molecule.definition = in_molecule.identity

                interactions_site = grna_tu_interactions.functionalComponents.create(site.displayId)

                interactions_site.name = site.name
                interactions_site.definition = site.identity

                interactions_repression = grna_tu_interactions.interactions.create('_represses_'.join([interactions_in_molecule.displayId, interactions_site.displayId]))
                interactions_repression.name = ' Represses '.join([in_molecule.name, site.name])
                interactions_repression.types = [SBO_INHIBITION]

                interactions_repressor = interactions_repression.participations.create(interactions_in_molecule.displayId)
                interactions_repressor.name = interactions_in_molecule.name
                interactions_repressor.roles = interactions_repressor.roles + [SBO_INHIBITOR]
                interactions_repressor.participant = interactions_in_molecule.identity

                interactions_repressed = interactions_repression.participations.create(interactions_site.displayId)
                interactions_repressed.name = interactions_site.name
                interactions_repressed.roles = interactions_repressed.roles + [SBO_INHIBITED]
                interactions_repressed.participant = interactions_site.identity

for network in networks:
    doc2.addModuleDefinition(network)
                        
doc2.write('generate_circuit_test2.xml')