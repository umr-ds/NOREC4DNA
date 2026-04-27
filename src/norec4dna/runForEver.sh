#!/bin/bash
echo $$ > runForEver.pid
exit_script() {
    echo "Ending Script."
    bash upload.sh
    rm runForEver.pid
    trap - SIGINT SIGTERM # clear the trap
    kill -- -$$ # Sends SIGTERM to child/sub processes
}
# react to kill -INT <PID>
trap exit_script SIGINT SIGTERM

echo "Running for ever:"
while true
do
    python3 -m norec4dna.simulator.simulator_new
done
