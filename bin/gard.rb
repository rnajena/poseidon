#!/usr/bin/env ruby

def adjust_gard_html(gard_html_file, html_kh_text, breakpoints_h, kh_insignificant_bp)
    gard_html_file_adjusted = File.open("#{gard_html_file.sub('.html','.adjusted.html')}",'w')
    f = File.open(gard_html_file,'r')
    f.each do |l|
      if l.start_with?('</DIV>')
        gard_html_file_adjusted << l.sub('</DIV></html></body></html>',"\n\n")
        # add additional content for KH test
        gard_html_file_adjusted << html_kh_text
        # add final breakpoints identified and if they were adjusted
        if kh_insignificant_bp
          gard_html_file_adjusted << '<p><u>Final breakpoints</u><br>(regardless of the KH-test (user defined) and adjusted to keep in-frame alignment, if necessary)</p>'
        else
          gard_html_file_adjusted << '<p><u>Final breakpoints</u><br>(with significant (adjp<0.01) KH-test and adjusted to keep in-frame alignment, if necessary)</p>'
        end
        if breakpoints_h.size > 0
          gard_html_file_adjusted << '<ul>'
          breakpoints_h.each do |bp_pos, p|
            gard_html_file_adjusted << "<li><b>#{bp_pos}</b></li>"
          end
          gard_html_file_adjusted << '</ul>'
        else
          gard_html_file_adjusted << '<p>No breakpoints with significant topological incongruence found.</p>'
        end
      else
        if l.start_with?('<!DOCTYPE')
          gard_html_file_adjusted << l.split('alignment <b>')[0] << 'alignment with' << l.split('</b> with')[1]
        else
          gard_html_file_adjusted << l
        end
      end
    end
    if kh_insignificant_bp
      gard_html_file_adjusted << "<p><b>ATTENTION</b>: as defined by the user break points are also used for further calculations even if the KH test showed no significant topological incongruence!</p>\n"
    end

    gard_html_file_adjusted << '</DIV></html></body></html>'

    f.close
    gard_html_file_adjusted.close

    #`cp #{f.path} #{f.path}.save` unless File.exists?("#{f.path}.original.save")
    #`mv #{gard_html_file_adjusted.path} #{f.path}`
    #`cp #{f.path} #{f.path}.adjusted.save` unless File.exists?("#{f.path}.adjusted.save")
end


def collect_params(output, html_text, kh_insignificant_bp)
    # collect: 1) breakpoints, 2) KH significance

    breakpoints = {}

    read = false

    file = File.open("gard_processor.log",'r')
    file.each do |line|
      if read

        if line.include?('|')
          html_text << '<tr><td>' << line.gsub('|','</td><td>') << '</td></tr>'
        else
          html_text << '</tbody></table>' if line.start_with?('At') && html_text.end_with?('</tr>')
          html_text << "\n<p>" << line.chomp << '</p>' if line.start_with?('At') || line.start_with?('Mean')
        end

        if line.include?('|')

          s = line.split('|')
          bp_pos = s[0].to_i
          lhs_adjp = s[2].to_f
          rhs_adjp = s[4].to_f

          # check if breakpoint is modulo 3, if not --> adjust
          bp_pos += 1 while bp_pos.modulo(3) != 0

          if lhs_adjp < 0.1 && rhs_adjp < 0.1
            # this breakpoint is significant!
            breakpoints[bp_pos] = 0.1 if lhs_adjp < 0.1 && rhs_adjp < 0.1
            breakpoints[bp_pos] = 0.05 if lhs_adjp < 0.05 && rhs_adjp < 0.05
            breakpoints[bp_pos] = 0.01 if lhs_adjp < 0.01 && rhs_adjp < 0.01
          else
            if kh_insignificant_bp
              breakpoints[bp_pos] = 1
            end
          end
        end
      end

      if line.start_with?('Breakpoint')
        read = true
        html_text << '<thead><tr><th>' << line.gsub('|','</th><th>') << '</th></tr></thead><tbody>'
      end
    end
    file.close

    # write out breakpoints to file for further usage
    bp = File.open('bp.tsv','w')
    breakpoints.each do |pos, significance|
      bp << "#{pos}\t#{significance}\n"
    end
    bp.close

    breakpoints
end

kh_insignificant_bp = ARGV[0].to_s.downcase == "true"
html_kh_text = '<p><u>KH-test</u></p><table>'
output = '.'

## COLLECT BREAKPOINT POSITIONS IF THERE ARE ANY
breakpoints = collect_params(output, html_kh_text, kh_insignificant_bp)
puts "We found #{breakpoints.size} significant break points for further analyses:\n\t\t#{breakpoints}\n"

gard_html_file = ARGV[1]

## ADJUST THE GARD OUTPUT HTML FILE
adjust_gard_html(gard_html_file, html_kh_text, breakpoints, kh_insignificant_bp)

