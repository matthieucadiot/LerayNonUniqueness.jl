#!/bin/bash
set -eu
exec "$1" "$2" worker "$3" "$4"
