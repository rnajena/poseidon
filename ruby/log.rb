#!/home/hoelzer/local/bin/ruby

timestamp = Time.now.to_i
`cp /mnt/fass2/poseidon-webserver-prost/poseidon_run.log /mnt/fass2/poseidon-webserver-prost/logs/#{timestamp}.poseidon_run.log`
`cp /mnt/fass2/poseidon-webserver-mahlzeit/poseidon_run.log /mnt/fass2/poseidon-webserver-mahlzeit/logs/#{timestamp}.poseidon_run.log`
`cp /mnt/fass2/poseidon-webserver-dessert/poseidon_run.log /mnt/fass2/poseidon-webserver-dessert/logs/#{timestamp}.poseidon_run.log`