mkfifo -m 600 logpipe
cat <> logpipe 1>&2 &
make create-capsule-sphere.tst
diff create-capsule-sphere/log create-capsule-sphere.ref

