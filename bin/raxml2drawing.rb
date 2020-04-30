#!/usr/bin/env ruby

output = File.open(ARGV[1], 'w')

File.open(ARGV[0], 'r').each do |newick|
  # get bootstrap values
  new_newick = newick
  newick.scan(/:[0-9]+\.[0-9]+\[[0-9]+\]/).each do |bootstrap|
    b = bootstrap.split('[')[1].sub(']','')
    n = bootstrap.split('[')[0]
    new_bootstrap = "#{b}#{n}"
    new_newick = new_newick.sub(bootstrap, new_bootstrap)
  end
  # round branch lengths
  newick.scan(/[0-9]+\.[0-9]+/).each do |length|
    f = length.to_f.round(4)
    new_newick = new_newick.sub(length, "#{f}")
  end
  output << new_newick
end
output.close