import os
import subprocess
import sys


def main() -> int:
    command = [
        sys.executable,
        "-m",
        "unittest",
        "discover",
        "-s",
        "tests/integration",
        "-v",
    ]
    env = os.environ.copy()
    env["RUN_LIVE_TESTS"] = "1"
    completed = subprocess.run(command, env=env)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
