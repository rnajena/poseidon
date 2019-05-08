#!/home/hoelzer/local/bin/ruby

# 1) summarize all CodeML files (ctl, mlc)
# 2) build also a table of the LRT tests
# 3) this need to be done for the full aln and also for fragments in case of fragmentation

require 'fileutils'

class CodemlHtml

  def initialize(type, codeml_html_out, html_index_file, tex_objects, frag_names)

    models = %w(M0 M1a M2a M7 M8 M8a)
    freqs = %w(F3X4 F1X4 F61)

    init_html(codeml_html_out, html_index_file, type, frag_names)

    if type == 'full_aln'
      codeml_html_out << "\n<h1>Results of CodeML analyses for the full alignment</h1>\n"
    else
      codeml_html_out << "\n<h1>Results of CodeML analyses for #{type.sub('_',' ')}</h1>\n"
    end

    # build table with all result files of the CodeML runs
    codeml_html_out << "<br><h2>CodeML parameter and result files</h2>\n"
    codeml_html_out << overview_codeml(models, freqs)

    # build a table of the LRT test parameters
    codeml_html_out << "<br><h2>Likelihood ratio tests between nested models</h2>\n"
    codeml_html_out << "<p>Green = LRT was significant at significance level 0.1</p><p>Red = LRT was not significant.</p>\n"
    codeml_html_out << lrt(tex_objects[type])


    codeml_html_out << '</div></div></body></html>'
    codeml_html_out.close

  end

  def lrt(tex_objects)

    html_table = "\n<table cellpadding=\"10\" border=\"1\"><thead><tr><th></th><th>M1a-M2a</th><th>M7-M8</th><th>M8a-M8</th></tr></thead><tbody>"

    tex_objects.each do |tex_object|
      freq = File.dirname(tex_object.output_file.path).reverse.split('/')[0].reverse
      html_table << "<tr><td><b>#{freq}</b></td>"

      lrt_m2a = tex_object.significance_m2a.lrt
      pvalue_m2a = tex_object.significance_m2a.pvalue.gsub('$','').sub('<','').to_f
      lnL_m1a = tex_object.significance_m2a.lnL0
      lnL_m2a = tex_object.significance_m2a.lnL1
      omega_percent_m2a = tex_object.significance_m2a.omega_percent
      omega_average_m2a = tex_object.significance_m2a.average_omega

      lrt_m8 = tex_object.significance_m8.lrt
      lrt_m8a = tex_object.significance_m8a.lrt
      pvalue_m8 = tex_object.significance_m8.pvalue.gsub('$','').sub('<','').to_f
      pvalue_m8a = tex_object.significance_m8a.pvalue.gsub('$','').sub('<','').to_f
      lnL_m7 = tex_object.significance_m8.lnL0
      lnL_m8a = tex_object.significance_m8a.lnL0
      lnL_m8 = tex_object.significance_m8.lnL1
      omega_percent_m8 = tex_object.significance_m8.omega_percent
      omega_average_m8 = tex_object.significance_m8.average_omega

      # write M1a-M2a parameters
      color = '#A9F5A9'
      comparison = '='
      if pvalue_m2a >= 0.1
        color = '#FA5858'
        omega_average_m2a = 'NA'
        omega_percent_m2a = 'NA'
      end
      comparison = '<' if pvalue_m2a == 0.001
      html_table << "<td bgcolor=\"#{color}\">lnL_<sub>M1a</sub> = #{lnL_m1a} <br> lnL_<sub>M2a</sub> = #{lnL_m2a} <br> LRT: &chi;<sup>2</sup> = #{lrt_m2a} <br> % sites with &omega; > 1 = #{omega_percent_m2a} <br> avg(&omega;) = #{omega_average_m2a} <br> p #{comparison} #{pvalue_m2a}</td>\n"

      # write M7-M8 parameters
      color = '#A9F5A9'
      comparison = '='
      if pvalue_m8 >= 0.1
        color = '#FA5858'
        omega_average_m8 = 'NA'
        omega_percent_m8 = 'NA'
      end
      comparison = '<' if pvalue_m8 == 0.001
      html_table << "<td bgcolor=\"#{color}\">lnL_<sub>M7</sub> = #{lnL_m7} <br> lnL_<sub>M8</sub> = #{lnL_m8} <br> LRT: &chi;<sup>2</sup> = #{lrt_m8} <br> % sites with &omega; > 1 = #{omega_percent_m8} <br> avg(&omega;) = #{omega_average_m8} <br> p #{comparison} #{pvalue_m8}</td>\n"

      # write M8a-M8 parameters
      color = '#A9F5A9'
      comparison = '='
      if pvalue_m8a >= 0.1
        color = '#FA5858'
        omega_average_m8 = 'NA'
        omega_percent_m8 = 'NA'
      end
      comparison = '<' if pvalue_m8a == 0.001
      html_table << "<td bgcolor=\"#{color}\">lnL_<sub>M8a</sub> = #{lnL_m8a} <br> lnL_<sub>M8</sub> = #{lnL_m8} <br> LRT: &chi;<sup>2</sup> = #{lrt_m8a} <br> % sites with &omega; > 1 = #{omega_percent_m8} <br> avg(&omega;) = #{omega_average_m8} <br> p #{comparison} #{pvalue_m8a}</td>\n"

      html_table << "</tr>\n"
    end

    html_table << "</tbody></table>\n"
    html_table
  end

  def overview_codeml(models, freqs)
    ctl_path = 'params/codeml_configs/'
    mlc_path = 'data/codeml/'

    overview_table = "\n<table cellpadding=\"10\" border=\"1\"><thead><tr><th></th>"
    models.each do |model|
      overview_table << "<th>#{model}</th>"
    end
    overview_table << "</tr></thead>\n<tbody>\n"
    freqs.each do |freq|
      overview_table << "<tr><td><b>#{freq}</b></td>"
      models.each do |model|
        overview_table << "<td><a href=\"#{mlc_path}/#{freq}/#{model}/codeml.mlc\">results</a> | <a href=\"#{ctl_path}/#{freq}/#{model}/#{freq}_#{model}_codeml.ctl\">params</a></td>\n"
      end
      overview_table << "</tr>\n"
    end
    overview_table << "</tbody></table>\n"
    overview_table
  end

  def init_html(html_out_file, html_index_file, type, frag_names)
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
      html_out_file << line.sub('<h4>full alignment</h4>',"<h4>#{type}</h4>").sub('145px','185px').sub('a:','a.').chomp << "\n" if write
      if line.start_with?('<h3>')
        # include link to the main output
        html_out_file << "<p><b><a class=\"link\" href=\"index.html\">Go back to Positive Selection Results</a></b></p>"
        if frag_names.size > 0
          html_out_file << "<details open=\"open\"><summary>Switch CodeML results</summary><ul>\n"
#          if type == 'full_aln'
#            html_out_file << "<li><a class=\"link\" href=\"codeml.html\">CodeML full aln</a></li>\n"
#          else
          unless type == 'full_aln'
            html_out_file << "<li><a class=\"link\" href=\"../full_aln/codeml.html\">CodeML full alignment</a></li>\n"
          end
          frag_names.each do |frag|
            if type != frag
              html_out_file << "<li><a class=\"link\" href=\"../#{frag}/codeml.html\">CodeML #{frag}</a></li>\n"
            end
          end
          html_out_file << "</ul></details>\n"
        end
      end
    end
    html_out_file << "</div></div><div id=\"maincontent\"><div id=\"innertube\"><br>\n"
    f.close
  end

end
