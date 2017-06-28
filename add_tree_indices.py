from pwgsresults.json_writer import JsonWriter as PwgsJson
from pwgsresults.result_loader import ResultLoader as PwgsResults
import json
import argparse

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('treesumm_fn')
  parser.add_argument('mutlist_fn')
  args = parser.parse_args()

  loader = PwgsResults(args.treesumm_fn, args.mutlist_fn, None)
  writer = PwgsJson(loader.dataset_name)
  writer.write_summaries(loader.tree_summary, loader.params, args.treesumm_fn)

main()
