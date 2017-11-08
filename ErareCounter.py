import molgenis, pprint
from configParser import ConfigParser

class ErareCounter:
    def __init__(self):
        config = ConfigParser().config
        self.genes = config['genes']
        self.labs = config['classification_columns']
        self.consensus_table = config['consensus_table']
        self.session = molgenis.Session(config['url'])
        self.session.login(config['account'], config['password'])

    def get_consensus_data(self):
        data = self.session.get(self.consensus_table, num=10000)
        labs = self.labs
        consensus = {}
        gene_stats = {}
        for variant in data:
            gene = variant['gene']['gene']
            if gene in self.genes:
                id = variant['id']
                cdna = variant['cDNA']
                try:
                    protein = variant['protein']
                except:
                    protein = ''
                try:
                    ref = variant['REF']
                except:
                    ref = ''
                try:
                    alt = variant['ALT']
                except:
                    alt = ''
                try:
                    chr = variant['#CHROM']
                except:
                    chr = ''
                try:
                    pos = variant['POS']
                except:
                    pos = ''
                if gene not in consensus:
                    gene_stats[gene] = {
                        'b': 0, 'p': 0, 'v': 0, 'no_consensus': 0, 'no_classifications': 0, 'one_classification': 0
                    }
                    consensus[gene] = [{'id': id,
                                        'cDNA': cdna,
                                        'chr': chr,
                                        'protein': protein,
                                        'ref': ref,
                                        'alt': alt,
                                        'pos': pos}
                                       ]
                else:
                    consensus[gene].append({'id': id,
                                            'cDNA': cdna,
                                            'chr': chr,
                                            'protein': protein,
                                            'ref': ref,
                                            'alt': alt,
                                            'pos': pos})
                labCount = 0
                classifications = []
                for lab in labs:
                    try:
                        classification = variant[lab]['classification']
                        labCount += 1
                        if classification == 'Benign':
                            classifications.append('b')
                        elif classification == 'Likely benign':
                            classifications.append('lb')
                        elif classification == "Pathogenic":
                            classifications.append('p')
                        elif classification == "Likely pathogenic":
                            classifications.append('lp')
                        else:
                            classifications.append('v')
                    except:
                        pass
                summary_counts = {'b': classifications.count('b') + classifications.count('lb'),
                                  'p': classifications.count('p') + classifications.count('lp'),
                                  'v': classifications.count('v')}
                complete_counts = {'b': classifications.count('b'), 'lb': classifications.count('lb'),
                                   'v': classifications.count('v'), 'lp': classifications.count('lp'),
                                   'p': classifications.count('p')}
                consensus[gene][-1]['counts'] = complete_counts

                consensus[gene][-1]['times_classified'] = labCount
                consensus[gene][-1]['status'] = 'Unsolved'
                # Variant should be classified at least twice in order to find consensus
                if labCount > 1:
                    lb_b = self.get_percentage(labCount, summary_counts['b'])
                    lp_p = self.get_percentage(labCount, summary_counts['p'])
                    b = self.get_percentage(labCount, complete_counts['b'])
                    lb = self.get_percentage(labCount, complete_counts['lb'])
                    p = self.get_percentage(labCount, complete_counts['p'])
                    lp = self.get_percentage(labCount, complete_counts['lp'])
                    v = self.get_percentage(labCount, summary_counts['v'])

                    # Variant is likely pathogenic when likely pathogenic is classified 50% or more
                    # and if the likely pathogenic classification is not equal to the vus classification
                    # and likely benign or benign is not classified
                    # OR
                    # Variant is classified pathogenic more than 50% and likely pathogenic and pathogenic are not equally
                    # classified and VUS is classified and likely benign or benign is not classified
                    if (lp >= 50 and lp != p and lp != v and lb_b == 0) or (p >= 50 and p != lp and v != 0 and lb_b == 0):
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Likely pathogenic'
                        consensus[gene][-1]['status'] = self.validate(lp)
                        gene_stats[gene]['p'] += 1

                    # Variant is classified pathogenic when pathogenic is classified 50% or more and the pathogenic
                    # classification is used more times than the likely pathogenic classification and the variant is never
                    # classified benign or likely benign (at this point, VUS can not be more than 0 since we checked for
                    # that in the previous check
                    elif p >= 50 and p != lp and lb_b == 0:
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Pathogenic'
                        consensus[gene][-1]['status'] = self.validate(p)
                        gene_stats[gene]['p'] += 1

                    # Variant is likely benign when likely benign is classified 50% or more
                    # and if the likely benign classification is not equal to the vus classification
                    # and likely pathogenic or pathogenic is not classified
                    # OR
                    # Variant is classified benign more than 50% and likely benign and benign are not equally
                    # classified and VUS is classified and likely pathogenic or pathogenic is not classified
                    elif (lb >= 50 and lb != b and lb != v and lp_p == 0) or (b >= 50 and b != lb and v != 0 and lp_p == 0):
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Likely benign'
                        consensus[gene][-1]['status'] = self.validate(lb)
                        gene_stats[gene]['b'] += 1

                    # Variant is classified benign when benign is classified 50% or more and the benign
                    # classification is used more times than the likely benign classification and the variant is never
                    # classified pathogenic or likely pathogenic (at this point, VUS can not be more than 0 since we checked
                    # for that in the previous check
                    elif b >= 50 and b != lb and lp_p == 0:
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Benign'
                        consensus[gene][-1]['status'] = self.validate(b)
                        gene_stats[gene]['b'] += 1

                    # Variant is classified likely benign when likely benign and benign together are classified more than
                    # 50% of the time. VUS should not be classified the same number of times as benign/likely benign and
                    # this variant is never classified as likely pathogenic and pathogenic
                    elif lb_b >= 50 and lb_b != v and lp_p == 0:
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Likely benign'
                        consensus[gene][-1]['status'] = self.validate(lb_b)
                        gene_stats[gene]['b'] += 1

                    # Variant is classified likely pathogenic when likely pathogenic and pathogenic together are classified
                    # more than 50% of the time. VUS should not be classified the same number of times as pathogenic/
                    # likely pathogenic and this variant is never classified as likely benign and benign
                    elif lp_p >= 50 and lp_p != v and lb_b == 0:
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Likely pathogenic'
                        consensus[gene][-1]['status'] = self.validate(lp_p)
                        gene_stats[gene]['p'] += 1

                    # Variant is classified VUS when VUS is classified 50% or more of the time. VUS should not be
                    # classified the same number of times as pathogenic/likely pathogenic and VUS should not be
                    # classified the same number of times as benign/likely benign
                    elif v >= 50 and v != lb_b and v != lp_p:
                        consensus[gene][-1]['consensus'] = 'Yes'
                        consensus[gene][-1]['classification'] = 'Uncertain significance (VOUS)'
                        consensus[gene][-1]['status'] = self.validate(v)
                        gene_stats[gene]['v'] += 1

                    # All other variants did not have consensus
                    else:
                        consensus[gene][-1]['classification'] = ''
                        gene_stats[gene]['no_consensus'] += 1
                        consensus[gene][-1]['consensus'] = 'No'

                # When a variant is classified less than 2 times, either the variant is never classified, or the variant
                # is classified by one lab
                else:
                    consensus[gene][-1]['classification'] = ''
                    if labCount == 0:
                        gene_stats[gene]['no_classifications'] += 1
                        consensus[gene][-1]['consensus'] = 'Not classified'
                    elif labCount == 1:
                        gene_stats[gene]['one_classification'] += 1
                        consensus[gene][-1]['consensus'] = 'No, classified by one lab'
            self.write_results(consensus, gene_stats)


    def validate(self, percentage):
        if percentage >= 75:
            return 'Validated'
        else:
            return 'Provisional'

    def get_percentage(self, length, count):
        return count * 100 / length

    def write_results(self, consensus, gene_stats):
        v_output = open('variant_output.csv', 'w')
        v_output.write(
            '"gene","id","chr","pos","ref","alt","cdna","protein","isConsensus",' +
            '"classification","status","times_classified","b","lb","v","lp","p"\n')
        for gene in consensus:
            for variant in consensus[gene]:
                v_output \
                    .write(
                    '"{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}","{}"\n'.format(gene,
                                                                                                               variant[
                                                                                                                   'id'],
                                                                                                               str(
                                                                                                                   variant[
                                                                                                                       'chr']),
                                                                                                               str(
                                                                                                                   variant[
                                                                                                                       'pos']),
                                                                                                               variant[
                                                                                                                   'ref'],
                                                                                                               variant[
                                                                                                                   'alt'],
                                                                                                               variant[
                                                                                                                   'cDNA'],
                                                                                                               variant[
                                                                                                                   'protein'],
                                                                                                               variant[
                                                                                                                   'consensus'],
                                                                                                               variant[
                                                                                                                   'classification'],
                                                                                                               variant[
                                                                                                                   'status'],
                                                                                                               variant[
                                                                                                                   'times_classified'],
                                                                                                               variant[
                                                                                                                   'counts'][
                                                                                                                   'b'],
                                                                                                               variant[
                                                                                                                   'counts'][
                                                                                                                   'lb'],
                                                                                                               variant[
                                                                                                                   'counts'][
                                                                                                                   'v'],
                                                                                                               variant[
                                                                                                                   'counts'][
                                                                                                                   'lp'],
                                                                                                               variant[
                                                                                                                   'counts'][
                                                                                                                   'p']))
        p_output = open('protein_output.csv', 'w')
        p_output.write('"GENE","BENIGN","PATHOGENIC","VOUS","NO CONSENSUS","NO CLASSIFICATIONS","ONE CLASSIFICATION"\n')
        for gene in gene_stats:
            p_output.write('"{}","{}","{}","{}","{}","{}","{}"\n'.format(gene, gene_stats[gene]['b'], gene_stats[gene]['p'],
                                                                    gene_stats[gene]['v'],
                                                                    gene_stats[gene]['no_consensus'],
                                                                    gene_stats[gene]['no_classifications'],
                                                                    gene_stats[gene]['one_classification']))
        v_output.close()
        p_output.close()


def main():
    erareCount = ErareCounter()
    erareCount.get_consensus_data()
    print('Done!')


if __name__ == '__main__':
    main()
