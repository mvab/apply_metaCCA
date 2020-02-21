#!/usr/bin/python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import requests
import os
import json
import pandas as pd


def parse_args():
    """Parse command line arguments"""
    parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('rsid_list', help="List of rsids to get proxies for")
    parser.add_argument('r_sq', nargs='?', default=0.8, help="Max r2 LD value for proxies")
    args = parser.parse_args()
    return args


def main():
    """ Main """
    args = parse_args()
    r_sq = float(args.r_sq)
    rsid_list = args.rsid_list

    if not os.path.isfile(rsid_list):
        print("Either file is missing or is not readable:", rsid_list)
    rsids = [rs.strip() for rs in open(rsid_list)]

    main_df = pd.DataFrame()
    for rsid in rsids:
        query_json = query_api(rsid)
        rsid_df = parse_json(query_json)
        main_df = pd.concat([main_df, rsid_df])

    main_df_subset = main_df[main_df['proxy_rsq'] > r_sq]
    main_df_subset.to_csv("ld_proxies_output.csv", index=False)

    print("# Summary:")
    print("# - Looked up {} rsids".format(len(set(rsids))))
    print("# - Returned {} proxies".format(main_df.shape[0]))
    print("# - Number of proxies that met {} LD r2 threshold: {}".format(r_sq, main_df_subset.shape[0]))
    print("# - Out of them {} map to unique original rsids ".format(len(set(main_df_subset['rsid_given']))))



def query_api(rsid):
    """Run proxy query for a specified rsid, return json """

    print(" > Looking up proxies for {}".format(rsid))
    headers = {'Content-Type': 'application/json'}
    data = json.dumps({"size": 10, "query": {"bool": {"filter": [{"term": {"target": "%s"}}]}}}) % rsid
    response = requests.get("http://140.238.83.192:9200/mrb-proxies/_search?pretty", headers=headers, data=data)

    return response.json()


def parse_json(query_json):
    """Parse json query output to extract individual proxies"""

    df = pd.DataFrame()
    hits = query_json["hits"]["hits"]

    if len(hits) != 0:
        print("Found {} rsid proxies".format(len(hits)))

        for i in list(range(0, len(hits))):
            # extract all rsid:proxy pairs for a given rsid
            hit = hits[i]
            rsid_given = hit['_source']['target']
            rsid_proxy = hit['_source']['proxy']
            proxy_rsq = hit['_source']['rsq']

            # add to pandas df
            temp = pd.DataFrame({'rsid_given': rsid_given,
                                 'rsid_proxy': rsid_proxy,
                                 'proxy_rsq': proxy_rsq}, index=[i])
            df = pd.concat([df, temp])
    else:
        print("===== No proxies found! =====")
    return df


if __name__ == '__main__':
    main()
