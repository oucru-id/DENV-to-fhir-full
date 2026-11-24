#!/usr/bin/env python3
"""Fetch an OAuth2 client-credentials token and save it outside the repository (~/.denvfhir)."""

import argparse
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

DEFAULT_HOME = Path(os.environ.get("DENVFHIR_HOME") or Path.home() / ".denvfhir")


def fetch_token(auth_base_url, client_id, client_secret):
    token_url = f"{auth_base_url.rstrip('/')}/protocol/openid-connect/token"
    form = urllib.parse.urlencode({
        "grant_type": "client_credentials",
        "client_id": client_id,
        "client_secret": client_secret,
        "scope": "openid",
    }).encode()
    req = urllib.request.Request(
        token_url, data=form,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        sys.exit(f"ERROR: token request failed (HTTP {exc.code}): {exc.read().decode('utf-8', 'replace')}")
    except urllib.error.URLError as exc:
        sys.exit(f"ERROR: cannot reach {token_url}: {exc.reason}")

    access_token = payload.get("access_token")
    if not access_token:
        sys.exit("ERROR: no access_token in response")
    now = datetime.now(timezone.utc)
    expires_in = int(payload.get("expires_in", 3600))
    return {
        "token_url": token_url,
        "issued_at": now.isoformat(),
        "expires_at": (now + timedelta(seconds=expires_in)).isoformat(),
        "access_token": access_token,
        "token_type": payload.get("token_type", "Bearer"),
        "expires_in": expires_in,
        "scope": payload.get("scope", ""),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auth_base_url", default=os.environ.get("DENVFHIR_AUTH_URL", ""))
    parser.add_argument("--client_id", default=os.environ.get("DENVFHIR_CLIENT_ID", ""))
    parser.add_argument("--client_secret", default=os.environ.get("DENVFHIR_CLIENT_SECRET", ""))
    parser.add_argument("--output", default=str(DEFAULT_HOME / "access_token.json"))
    args = parser.parse_args()

    if not (args.auth_base_url and args.client_id and args.client_secret):
        sys.exit("ERROR: provide --auth_base_url, --client_id, and --client_secret (or set DENVFHIR_* env vars)")

    token_data = fetch_token(args.auth_base_url, args.client_id, args.client_secret)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(token_data, indent=2), encoding="utf-8")
    os.chmod(out, 0o600)
    print(f"Token saved to {out} (expires {token_data['expires_at']})")


if __name__ == "__main__":
    main()
