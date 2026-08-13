"""
Minimal connectivity check for the remote Fluent gRPC session running on
the Windows PC. Does nothing else -- confirms the network path and
authentication work before anything in fluent_solve.py / fluent_batch.py
is attempted.

Usage (pick ONE of the two connection methods -- see FLUENT_REMOTE.md for
how to obtain these values from the Windows machine):

    # via a server-info file copied from the Windows machine
    python3 test_connection.py --server-info-file path/to/server_info.txt

    # via explicit IP/port/password
    python3 test_connection.py --ip 192.168.1.50 --port 12345 --password abc123

VERIFIED vs NOT: connect_to_fluent's signature below (ip/port/address/
server_info_file_name/password/allow_remote_host/...) was checked directly
against the installed ansys-fluent-core package source in this environment
(pip install ansys-fluent-core), so the argument names are real, not
guessed. What is NOT verified is an actual live connection -- that can
only happen once Fluent's gRPC server is reachable from this machine,
which is exactly what this script exists to confirm, on your machine.
"""

from __future__ import annotations

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Confirm connectivity to a remote Fluent gRPC session.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--server-info-file", help="Path to the server-info file written by Fluent on the Windows PC.")
    group.add_argument("--ip", help="IP address of the Windows PC running Fluent.")
    parser.add_argument("--port", type=int, help="gRPC port (required with --ip).")
    parser.add_argument("--password", help="gRPC password (required with --ip; also usable with --server-info-file "
                                            "if the file doesn't embed one).")
    args = parser.parse_args()

    if args.ip and args.port is None:
        parser.error("--port is required when using --ip")

    import ansys.fluent.core as pyfluent

    print("Connecting to Fluent...")
    try:
        if args.server_info_file:
            session = pyfluent.connect_to_fluent(
                server_info_file_name=args.server_info_file,
                password=args.password,
                allow_remote_host=True,
            )
        else:
            session = pyfluent.connect_to_fluent(
                ip=args.ip,
                port=args.port,
                password=args.password,
                allow_remote_host=True,
            )
    except Exception as e:
        print(f"FAILED to connect: {type(e).__name__}: {e}")
        print()
        print("Common causes: Fluent's gRPC server isn't started/exposed on the Windows PC yet, "
              "a firewall is blocking the port, the IP/port/password don't match what Fluent printed, "
              "or the server-info file is stale (each Fluent launch gets a new port/password). "
              "See FLUENT_REMOTE.md.")
        sys.exit(1)

    try:
        version = session.get_fluent_version()
        print(f"Connected. Fluent version reported: {version}")
        print(f"Session active: {session.is_active()}")
    finally:
        session.exit()
        print("Session closed cleanly.")


if __name__ == "__main__":
    main()
