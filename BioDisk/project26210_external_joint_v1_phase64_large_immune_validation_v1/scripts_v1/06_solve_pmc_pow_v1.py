#!/usr/bin/env python3
import hashlib
import re
import sys


def main(html_path):
    text = open(html_path, encoding="utf-8").read()
    challenge_match = re.search(r'const POW_CHALLENGE = "([^"]+)"', text)
    difficulty_match = re.search(r'const POW_DIFFICULTY = "(\d+)"', text)
    cookie_name_match = re.search(r'const POW_COOKIE_NAME = "([^"]+)"', text)
    if not challenge_match or not difficulty_match or not cookie_name_match:
        raise SystemExit("PMC PoW components were not found")
    challenge = challenge_match.group(1)
    difficulty = int(difficulty_match.group(1))
    prefix = "0" * difficulty
    nonce = 0
    while True:
        digest = hashlib.sha256(f"{challenge}{nonce}".encode()).hexdigest()
        if digest.startswith(prefix):
            print(f"{cookie_name_match.group(1)}={challenge},{nonce}")
            print(f"hash={digest}", file=sys.stderr)
            return
        nonce += 1


if __name__ == "__main__":
    main(sys.argv[1])

