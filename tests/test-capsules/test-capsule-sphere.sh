mkfifo -m 600 /tmp/logpipe
unbuffer cat /tmp/logpipe &
make create-capsule-sphere.tst
diff create-capsule-sphere/log create-capsule-sphere.ref

