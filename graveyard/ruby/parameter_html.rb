#!/usr/bin/env ruby

# 1) summarize all parameter files

class ParameterHtml

  PARAM_TXT = '../tools/tools.txt'

  def initialize(parameter_html_out, html_index_file, frag_names, timestamp, version, parameter_strings)

    param_files = {}
    parameter_strings.each do |type, parameter_text|
      out = File.open("#{File.dirname(html_index_file).sub('full_aln',type)}/params/#{type}_params.txt",'w')
      out << parameter_text
      out.close
      param_files[type] = out.path
    end

    init_html(parameter_html_out, html_index_file)

    parameter_html_out << "\n<h1>PoSeiDon (v#{version}) parameter settings and executed commands</h1>\n"

    zip = "http://www.rna.uni-jena.de/poseidon-results/#{timestamp}/#{timestamp}.zip"
    parameter_html_out << "\n<p><b><a href=\"#{zip}\">Download all results</a></b></p>"

    parameter_html_out << "\n<br><h2>Commands</h2>"
    parameter_html_out << parameter_files(parameter_strings, param_files)

    param_file = File.open(PARAM_TXT,'r')
    parameter_html_out << "\n<br><h2>Tools executed by PoSeiDon</h2>\n<ul>"
    param_file.each do |param_line|
      s = param_line.split(',')
      tool = s[0]
      version = "(#{s[1]})"
      version = '' if version == '(.)'
      author = ", #{s[2]}"
      author = '' if author == ', .'
      pmid = s[3].chomp
      pmid = '' if pmid == '.'
      if pmid.length > 0
        pmid = "<a target=\"_blank\" href=\"https://www.ncbi.nlm.nih.gov/pubmed/#{pmid}\">#{pmid}</a>"
      end
      parameter_html_out << "<li>#{tool} #{version}#{author} #{pmid}</li>"
    end
    parameter_html_out << "</ul>\n"
    param_file.close

    parameter_html_out << "\n<br><h2>The pipeline</h2>\n"
    parameter_html_out << "<a href=\"../src/pipeline_landscape.pdf\"><img width=\"1200px\" src=\"../src/pipeline_landscape.png\"></a>\n"


    parameter_html_out << '</div></div></body></html>'
    parameter_html_out.close

  end

  def parameter_files(parameter_strings, param_files)

    html = "<ul>\n"
    param_files.each do |type, file_path|
      html << "<li><a href=\"../#{type}/params/#{type}_params.txt\">Download commands #{type}</a></li>\n"
    end
    html << "</ul>\n\n"

    html << "<p><a href=\"codeml.html\">Get configuration files used to run CodeML</a></p>\n"

    parameter_strings.each do |type, text|
      html << "<details><summary>\nCommands #{type}\n</summary>\n<code></br>#{text.gsub("\n",'</br>')}</br></code></details>\n"
    end
    html
  end

  def init_html(html_out_file, html_index_file)
    f = File.open(html_index_file,'r')
    write = true
    refactor_ul = false
    f.each do |line|
      refactor_ul = true if line.include?('ul {')
      write = false if line.include?('<table><tr>')
      if line.include?('html body{')
        html_out_file << "\tp.alert {
    width: 800px;
    margin: 2px;
    #border-left: 6px solid #0044cc;
    background-color: #ccffff;
\t}\n\n
    details {
      cursor: pointer;
    }\n
    code {
    background: hsl(220, 80%, 90%);
}\n\n"
      end
      if refactor_ul && line.include?('padding:')
        html_out_file << line.sub('padding: 0','padding: 4')
        refactor_ul = false if line.include?('}')
        next
      end
      html_out_file << '#' if line.include?('overflow: scroll; /*Disable scrollbars')
      html_out_file << line.sub('<h4>full alignment</h4>','<h4>Parameters</h4>').sub('145px','185px').sub('a:','a.').chomp << "\n" if write
      if line.start_with?('<h3>')
        # include link to the main output
        html_out_file << "<p><b><a class=\"link\" href=\"index.html\">Go back to Positive Selection Results</a></b></p>"
      end
    end
    html_out_file << "</div></div><div id=\"maincontent\"><div id=\"innertube\"><br>\n"
    f.close
  end


end
