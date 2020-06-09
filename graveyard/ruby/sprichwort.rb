#!/usr/bin/env ruby

require 'open-uri'

class String
  # colorization
  def colorize(color_code)
    "\e[#{color_code}m#{self}\e[0m"
  end

  def red
    colorize(31)
  end

  def green
    colorize(32)
  end

  def yellow
    colorize(33)
  end

  def blue
    colorize(34)
  end

  def pink
    colorize(35)
  end

  def light_blue
    colorize(36)
  end
end

#url = %w(http://sprichwort.gener.at/or/ http://proverb.gener.at/or/).sample
url = 'http://sprichwort.gener.at/or/'

puts "\nTip of the Day:\t".yellow + "#{open(url).read.scan(/spwort\">.*</)[0].sub('spwort">','').sub('<','')}\n\n".light_blue





