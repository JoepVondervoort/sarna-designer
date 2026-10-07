"""Launch the local workbench: python -m sarna.app --port 8001."""

import argparse
from pathlib import Path

from sarna.app.service import DesignService
from sarna.app.web import AppRouter, make_server


def main() -> None:
    """Launch a loopback-only app with runs in a separate per-user data directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8001)
    parser.add_argument("--data-dir", type=Path, default=Path.home() / ".sarna-designer")
    args = parser.parse_args()
    if not 1 <= args.port <= 65535:
        parser.error("port must be between 1 and 65535")
    router = AppRouter(DesignService(args.data_dir))
    try:
        server = make_server(router, args.port)
    except OSError as error:
        router.close()
        parser.exit(1, f"Could not start the local app: {error}\nTry another --port.\n")
    print(f"Open http://localhost:{args.port}/", flush=True)
    print(f"Local runs: {args.data_dir / 'runs.sqlite3'}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        router.close()


if __name__ == "__main__":
    main()
