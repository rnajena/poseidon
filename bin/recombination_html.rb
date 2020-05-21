#!/usr/bin/env ruby

# 1) if recombination occured, show summary of GARD results
# 2) make user aware of recombination in his alignment
# 3) show the different tree topologies next to each other, for the moment just show the NT trees

require 'fileutils'

class RecombinationHtml

  def initialize(html_recomb_out, gard_html_file, full_tree_nt, fragment_trees_nt, html_index_out)

    html_out_file = File.open(html_recomb_out,'w')

    init_html(html_out_file, html_index_out)

    html_out_file << "\n<h1>Attention! PoSeiDon detected recombination in your alignment</h1><p class=\"alert\">Before you proceed, please be aware that
recombination events can have a profound impact on the calculated tree topologies and the positive selection detection. Especially, all results based on
the analyses of short fragments resulting from putative recombination events should be further investigated with suspicion.</p>\n\n"

    @use_insignificant = false
    gard_plot_html = select_gard_plot_html(gard_html_file)
    gard_kh_html = select_gard_kh_html(gard_html_file)

    if @use_insignificant
      html_out_file << "<br><p class=\"alert\">Please be aware that you defined to use break points for further calculations even if the KH test showed no significant topological incongruence! (\"Use insignificant break points\")</p>"
    end

    html_out_file << "<br><p class=\"alert\">Please use the navigation or click on the alignment links below to <a href=\"index.html\">check out the whole positive selection analysis</a>.</p>\n\n"

    html_out_file << "<br><h2>Recombination test with GARD</h2>\n\n"
    html_out_file << "<table cellpadding=\"20\" border=\"1\"><tr>\n"
    html_out_file << '<td>' << gard_plot_html << "<a href=\"data/gard.adjusted.html\">Detailed recombination output</a>" << '</td>'
    html_out_file << '<td>' << gard_kh_html << '</td>'
    html_out_file << "</tr></table>\n\n"

    html_out_file << "<br><h2>Resulting tree topologies</h2>\n<p>Trees on nucleotide level based on the full alignment and alignment fragments are shown.</p>\n"

    tree_counter = 0
    html_out_file << "<table cellpadding=\"20\" border=\"1\"><tr>\n"
#    html_out_file << "<td><b>Full alignment</b><br><br><a target=\"_blank\" href=\"data/#{File.basename(full_tree_nt).sub('.corrected.','.')}.scale.svg\">
    html_out_file << "<td><b>Full alignment</b><br><br><a target=\"_blank\" href=\"data/#{File.basename(full_tree_nt)}.corrected.scale.svg\">
<img src=\"data/#{File.basename(full_tree_nt)}.corrected.png\"></a><br>
<a href=\"index.html\">Go to alignment</a></td>"
    tree_counter += 1

    fragment_trees_nt.each do |fragment, tree|
      puts tree
      header = "Fragment #{fragment.split('_')[1]}"
      tree_counter += 1
#      html_out_file << '<td>' << "<b>#{header}</b><br><br><a target=\"_blank\" href=\"../#{fragment}/data/#{File.basename(tree).sub('.corrected.','.')}.scale.svg\">
      html_out_file << '<td>' << "<b>#{header}</b><br><br><a target=\"_blank\" href=\"../#{fragment}/data/#{File.basename(tree)}.corrected.scale.svg\">
<img src=\"../#{fragment}/data/#{File.basename(tree)}.corrected.png\"></a><br>
<a href=\"../#{fragment}/index.html\">Go to alignment</a></td>"
      if tree_counter == 2
        tree_counter = 0
        html_out_file << '</tr><tr>'
      end
    end
    html_out_file << "</tr></table>\n\n"

    html_out_file << '</div></div></body></html>'
    html_out_file.close

  end

  def init_html(html_out_file, html_index_file)
    f = File.open(html_index_file,'r')
    write = true
    f.each do |line|
      write = false if line.include?('<table><tr>')
      if line.include?('html body{')
        html_out_file << "\tp.alert {
    width: 800px;
    margin: 2px;
		padding: 2px;
    border-left: 6px solid #0044cc;
    #background-color: #ccffff;
\t}\n\n"
      end
      html_out_file << '#' if line.include?('overflow: scroll; /*Disable scrollbars')
      html_out_file << line.split('<h4>')[0].sub('145px','185px').sub('a:','a.').chomp << "\n" if write
      if line.start_with?('<h3>')
        # include link to the main output
        html_out_file << "<b><a class=\"link\" href=\"index.html\">Go to Positive Selection Results</a></b>"
      end
    end
    html_out_file << "</div></div><div id=\"maincontent\"><div id=\"innertube\"><br>\n"
    f.close
  end

  def select_gard_plot_html(file)
    gard_plot_html = ''
    read_plot = true
    f = File.open(file,'r')
    f.each do |l|
      read_plot = false if l.start_with?('</DIV><DIV CLASS=')
      if read_plot
        if l.include?('<u>Results</u><p>')
          gard_plot_html << l.split('<u>Results</u><p>')[1]
        else
          gard_plot_html << l
        end
      end
    end
    f.close
    gard_plot_html
  end

  def select_gard_kh_html(file)
    gard_kh_html = ''
    read_kh = false
    f = File.open(file,'r')
    f.each do |l|
      read_kh = true if l.start_with?('<p><u>KH-test')
      if read_kh
        @use_insignificant = true if l.include?('ATTENTION')
        gard_kh_html << l.split('<p><u>Final')[0] unless l.start_with?('</DIV>')
      end
    end
    f.close
    gard_kh_html
  end

end

#${gard_html} ${full_nt_tree} '${frag_nt_trees}'

html_index_out = ARGV[0]
gard_html_file = ARGV[1]
full_tree_nt = ARGV[2]

#[/home/martin/git/poseidon/work/b5/67a3f6987dce7bd1ad3efac06d8f14/bats_mx1_fragment_2_nt.raxml, /home/martin/git/poseidon/work/00/3c80a37ae31d549ad32e6eb4f6ac29/bats_mx1_fragment_3_nt.raxml, /home/martin/git/poseidon/work/c4/9fb484c6f2af416a59ea91043dbf3f/bats_mx1_fragment_1_nt.raxml]
fragment_trees_nt = {}
num_frags = ARGV[3].split(' ').size
num_frags.times.each do |i|
  ARGV[3].split(' ').each do |frag_tree|
    frag_id = 'fragment_'+frag_tree.sub('_nt.raxml','').split('fragment_')[1]
    fragment_trees_nt[frag_id] = frag_tree if "fragment_#{i+1}" == frag_id
  end
end

RecombinationHtml.new('html/full_aln/recomb.html', gard_html_file, full_tree_nt, fragment_trees_nt, html_index_out)