require 'logger'
require 'awesome_print'
require 'tempfile'
require 'bio'
require 'trollop'

$log = Logger.new STDERR

$opts = Trollop::options do
  opt :vcf, "VCF file to process", :type => :string, :required => false, :default => nil
  opt :fasta, "Reference fasta file", :type => :string, :required => false, :default => nil
  opt :method, "Method (fimo or tfscan)", :type => :string, :required => false, :default => "tfscan"
  opt :db, "DB to use (fimo only)", :type => :string, :required => false, :default => "data/HOCOMOCOv9_AD_MEME.txt"
  opt :window, "Window for transcription factor prediction", 
               :type => :int, 
               :required => false, :default => 15

  opt :reference, "File with reference sequence", :type => :string, :required => false, :default => nil
  opt :mutant, "File with mutant sequence", :type => :string, :required => false, :default => nil
end

def fimo(seq, database)
  $log.info "FIMO for: #{seq}"

  tmp = Tempfile.new 'fimo'
  tmp.puts ">seq\n#{seq}"
  tmp.close

  r = []

  IO.popen("fimo --norc --text --bgfile motif-file #{database} #{tmp.path} 2> /dev/null") do |io|
    io.each_line do |line|
      next if line[0] == "#"
      pattern, start, stop, pvalue = line.split.values_at(0, 2, 3, 6)
      start = start.to_i
      stop = stop.to_i
      next if start > $opts[:window] or stop < $opts[:window]

      r << [pattern, pvalue]
    end

    io.close
    throw("Failed to run fimo!") if $?.to_i != 0
  end

  return r
end

def tfscan(seq, database)
  $log.info "TFSCAN for: #{seq}"

  tmp = Tempfile.new 'tfscan'
  tmp.puts seq
  tmp.close

  r = []

  IO.popen("tfscan -sequence #{tmp.path} -menu V -mismatch 0 -filter 2> /dev/null") do |io|
    io.each_line do |line|
      next if line[0] == "#"
      s = line.split
      next if s.size < 6
      start = s[0].to_i
      stop = s[1].to_i
      pattern = s[4..-2] * " "
      next if start > $opts[:window] or stop < $opts[:window]

      r << [pattern, nil]
    end

    io.close
    throw("Failed to run fimo!") if $?.to_i != 0
  end

  r
end

def diff_pvalue(a, b)
  x = a.map(&:first) - b.map(&:first)
  a.find_all { |u, v| x.include?(u) }
end

def run_predictions(algo, seq, alt, db)
  ref = send(algo.to_s, seq.upcase, db)
  variant = send(algo.to_s, alt.upcase, db)

  gain = diff_pvalue(variant, ref)
  loss = diff_pvalue(ref, variant)

  [gain, loss]
end

def seq_with_mutation(s, a, p)
  front = s.slice(0, p)
  back = s[(p+1)..-1]

  front + "[" + s[p] + "/" + a + "]" + back
end

if $opts[:vcf] != nil
  throw("Need a reference file!") if $opts[:fasta] == nil
  File.foreach($opts[:vcf]) do |line|
    next if line[0] == "#"
    chr, pos, name, ref, alt = line.split("\t").values_at(0, 1, 2, 3, 4)
    throw("Not a SNP!") if ref.length != 1
    $log.info "SNP = #{name}"
    pos = pos.to_i
    alt = alt.split ","
    fa = %x{samtools faidx #{$opts[:fasta]} #{chr}:#{pos - $opts[:window]}-#{pos + $opts[:window]}}
    fa = Bio::FastaFormat.new(fa).seq
    fa = Bio::Sequence::NA.new(fa)
    rev = fa.reverse_complement.seq.upcase
    fa = fa.seq.upcase
    alt.each do |a|
      throw("Not a SNP!") if a.length != 1
      $log.info "Predictions for SNP in #{chr}:#{pos}"
      $log.info "Sequence is: #{seq_with_mutation(fa, a, $opts[:window]).upcase}"

      a_forward = fa.dup
      a_forward[$opts[:window]] = a
      gain, loss = run_predictions($opts[:method], fa, a_forward, $opts[:db])

      gain.each { |u| puts [chr, pos, name, ref, a, "+", "gain", u.first] * "\t" }
      loss.each { |u| puts [chr, pos, name, ref, a, "+", "loss", u.first] * "\t" }

      a_revcomp = Bio::Sequence::NA.new(a_forward).reverse_complement.seq
      $log.info "Sequence is: #{seq_with_mutation(rev, a_revcomp[$opts[:window]], $opts[:window]).upcase}"
      gain, loss = run_predictions($opts[:method], rev, a_revcomp, $opts[:db])

      gain.each { |u| puts [chr, pos, name, ref, a, "-", "gain", u].flatten * "\t" }
      loss.each { |u| puts [chr, pos, name, ref, a, "-", "loss", u].flatten * "\t" }
    end
  end
elsif $opts[:reference] != nil
  gain, loss = run_predictions($opts[:method], $opts[:reference], $opts[:mutant], $opts[:db])

  gain.each { |u| puts ["", "", "", "", "", "", "+", "gain", u.first] * "\t" }
  loss.each { |u| puts ["", "", "", "", "", "", "+", "loss", u.first] * "\t" }

  reference_revcomp = Bio::Sequence::NA.new($opts[:reference]).reverse_complement.seq
  mutant_revcomp = Bio::Sequence::NA.new($opts[:mutant]).reverse_complement.seq

  gain, loss = run_predictions($opts[:method], reference_revcomp, mutant_revcomp, $opts[:db])

  gain.each { |u| puts ["", "", "", "", "", "", "-", "gain", u].flatten * "\t" }
  loss.each { |u| puts ["", "", "", "", "", "", "-", "loss", u].flatten * "\t" }
else
  throw "Need either VCF with reference genome, or sequence and mutant inself!"
end
