"""
Functions parsing blast results
"""


def survey_blast_result(result, from_trimmed=True):
    import json

    def get_query_title(query):
        return query['report']['results']['search']['query_title'].split(' ')[0]
    if isinstance(result, str):
        with open(result) as handle:
            result = json.load(handle)
    if from_trimmed:
        return result

    if 'BlastOutput2' in result.keys():
        result = result['BlastOutput2']
    res = {}
    for query in result:
        query_title = get_query_title(query)
        res[query_title] = {'hit_num': len(query['report']['results']['search']['hits']),
                            'query_len': query['report']['results']['search']['query_len']}
        if res[query_title]['hit_num'] > 0:
            res[query_title]['hits'] = [{
                'accession': hit['description'][0]['accession'],
                'title': hit['description'][0]['title'],
                'best_hsp': hit['hsps'][0],
                'best_hsp_coverage': (abs(hit['hsps'][0]['query_to'] - hit['hsps'][0]['query_from']) + 1)/res[query_title]['query_len'],
                'best_hsp_align_len': hit['hsps'][0]['align_len'],
                'best_hsp_pct_identity': hit['hsps'][0]['identity']/hit['hsps'][0]['align_len']
            } for ix,hit in enumerate(query['report']['results']['search']['hits']) if ix <= 3]
    return res


def trim_blastn_results(results):
    """Due to large original blastn results, discard unnecessary entries"""
    res = survey_blast_result(results, from_trimmed=False)
    from pathlib import Path
    output_path = str(Path(results).parent) + '/' + str(Path(results).stem) + '_trimmed.json'
    import json
    with open(output_path, 'w') as handle:
        json.dump(obj=res, fp=handle)


def trim_blastn_results_main(data_root=None):
    from os import environ
    if data_root is None:
        if 'DATA_PATH' not in environ:
            raise EnvironmentError('Please indicate the root to data folder')
        else:
            root = environ['DATA_PATH'] + 'bglmm/'
    else:
        root = data_root + 'bglmm/'

    blast_res = ['fj_silva_aligned_16s.json', 'fj_silva_aligned_human.json',
                 'fj_silva_failed_16s.json', 'fj_silva_failed_human.json']

    for res in blast_res:
        trim_blastn_results(root + res)


if __name__ == '__main__':
    trim_blastn_results_main()
