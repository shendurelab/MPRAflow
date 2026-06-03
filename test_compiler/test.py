import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

def main():
    import mpranalyze_compiler as c
    print(c.parse("minimal_test_input.tsv"))

if __name__=="__main__":
    main()
