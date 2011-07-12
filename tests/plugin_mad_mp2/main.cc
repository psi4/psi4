<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
  
  

  


  

  <head>
    <title>
      /trunk/tests/plugin_mad_mp2/main.cc – Psi4
    </title>
        <link rel="search" href="/trac/search" />
        <link rel="help" href="/trac/wiki/TracGuide" />
        <link rel="alternate" href="/trac/browser/trunk/tests/plugin_mad_mp2/main.cc?format=txt" type="text/plain" title="Plain Text" /><link rel="alternate" href="/trac/export/1983/trunk/tests/plugin_mad_mp2/main.cc" type="text/x-c++src; charset=iso-8859-15" title="Original Format" />
        <link rel="up" href="/trac/browser/trunk/tests/plugin_mad_mp2" title="Parent directory" />
        <link rel="start" href="/trac/wiki" />
        <link rel="stylesheet" href="/trac/chrome/common/css/trac.css" type="text/css" /><link rel="stylesheet" href="/trac/chrome/common/css/code.css" type="text/css" /><link rel="stylesheet" href="/trac/pygments/trac.css" type="text/css" /><link rel="stylesheet" href="/trac/chrome/common/css/browser.css" type="text/css" /><link rel="stylesheet" href="/trac/chrome/googlemap/tracgooglemap.css" type="text/css" />
        <link rel="google-key" href="" class="google-key" title="ABQIAAAAEm-sB5Hq7DJNrZXGU2EXlxTthLP7MnSeoBrJ-Mt6ZcDfmRQM3RTqr4vj7oIAm9RhoGjV8d0hktaebQ" />
        <link rel="shortcut icon" href="/trac/chrome/site/logo.png" type="image/png" />
        <link rel="icon" href="/trac/chrome/site/logo.png" type="image/png" />
      <link type="application/opensearchdescription+xml" rel="search" href="/trac/search/opensearch" title="Search Psi4" />
    <script type="text/javascript" src="/trac/chrome/common/js/jquery.js"></script><script type="text/javascript" src="/trac/chrome/common/js/trac.js"></script><script type="text/javascript" src="/trac/chrome/common/js/search.js"></script><script type="text/javascript" src="/trac/chrome/googlemap/tracgooglemap.js"></script>
    <!--[if lt IE 7]>
    <script type="text/javascript" src="/trac/chrome/common/js/ie_pre7_hacks.js"></script>
    <![endif]-->
    <script type="text/javascript">
      jQuery(document).ready(function($) {
        $(".trac-toggledeleted").show().click(function() {
                  $(this).siblings().find(".trac-deleted").toggle();
                  return false;
        }).click();
        $("#jumploc input").hide();
        $("#jumploc select").change(function () {
          this.parentNode.parentNode.submit();
        });
      });
    </script>
  </head>
  <body>
    <div id="banner">
      <div id="header">
        <a id="logo" href="/trac/wiki"><img src="/trac/chrome/site/PSI4_3.png" alt="Welcome to the Psi4 Trac page" height="120" width="480" /></a>
      </div>
      <form id="search" action="/trac/search" method="get">
        <div>
          <label for="proj-search">Search:</label>
          <input type="text" id="proj-search" name="q" size="18" value="" />
          <input type="submit" value="Search" />
        </div>
      </form>
      <div id="metanav" class="nav">
    <ul>
      <li class="first">logged in as rparrish6</li><li><a href="/trac/logout">Logout</a></li><li><a href="/trac/prefs">Preferences</a></li><li><a href="/trac/wiki/TracGuide">Help/Guide</a></li><li class="last"><a href="/trac/about">About Trac</a></li>
    </ul>
  </div>
    </div>
    <div id="mainnav" class="nav">
    <ul>
      <li class="first"><a href="/trac/wiki">Wiki</a></li><li><a href="/trac/timeline">Timeline</a></li><li><a href="/trac/roadmap">Roadmap</a></li><li class="active"><a href="/trac/log/trunk?rev=latest">Browse Source</a></li><li><a href="/trac/report">View Tickets</a></li><li><a href="/trac/newticket">New Ticket</a></li><li><a href="/trac/search">Search</a></li><li><a href="/trac/doxygen">Doxygen</a></li><li><a href="/trac/admin" title="Administration">Admin</a></li><li class="last"><a href="/trac/discussion">Forum</a></li>
    </ul>
  </div>
    <div id="main">
      <div id="ctxtnav" class="nav">
        <h2>Context Navigation</h2>
          <ul>
              <li class="first"><a href="/trac/changeset/1970/trunk/tests/plugin_mad_mp2/main.cc">Last Change</a></li><li><a href="/trac/browser/trunk/tests/plugin_mad_mp2/main.cc?annotate=blame&amp;rev=1970" title="Annotate each line with the last changed revision (this can be time consuming...)">Annotate</a></li><li class="last"><a href="/trac/log/trunk/tests/plugin_mad_mp2/main.cc">Revision Log</a></li>
          </ul>
        <hr />
      </div>
    <div id="content" class="browser">
      <h1>
    <a class="pathentry first" title="Go to root directory" href="/trac/browser">root</a><span class="pathentry sep">/</span><a class="pathentry" title="View trunk" href="/trac/browser/trunk">trunk</a><span class="pathentry sep">/</span><a class="pathentry" title="View tests" href="/trac/browser/trunk/tests">tests</a><span class="pathentry sep">/</span><a class="pathentry" title="View plugin_mad_mp2" href="/trac/browser/trunk/tests/plugin_mad_mp2">plugin_mad_mp2</a><span class="pathentry sep">/</span><a class="pathentry" title="View main.cc" href="/trac/browser/trunk/tests/plugin_mad_mp2/main.cc">main.cc</a>
    <br style="clear: both" />
  </h1>
      <div id="jumprev">
        <form action="" method="get">
          <div>
            <label for="rev">
              View revision:</label>
            <input type="text" id="rev" name="rev" size="6" />
          </div>
        </form>
      </div>
      <div id="jumploc">
        <form action="" method="get">
          <div class="buttons">
            <label for="preselected">Visit:</label>
            <select id="preselected" name="preselected">
              <option selected="selected"></option>
              <optgroup label="branches">
                <option value="/trac/browser/trunk">trunk</option><option value="/trac/browser/branches/madpsi4">branches/madpsi4</option><option value="/trac/browser/branches/psi4-alpha-0">branches/psi4-alpha-0</option>
              </optgroup>
            </select>
            <input type="submit" value="Go!" title="Jump to the chosen preselected path" />
          </div>
        </form>
      </div>
      <table id="info" summary="Revision info">
        <tr>
          <th scope="col">
            Revision <a href="/trac/changeset/1970">1970</a>, <span title="1622 bytes">1.6 KB</span>
            (checked in by jturney, <a class="timeline" href="/trac/timeline?from=2011-07-09T20%3A25%3A37-0400&amp;precision=second" title="2011-07-09T20:25:37-0400 in Timeline">3 days</a> ago)
          </th>
        </tr>
        <tr>
          <td class="message searchable">
              <p>
Rob's new dfmp2 plugin to be made MAD(NESS)<br />
</p>
          </td>
        </tr>
      </table>
      <div id="preview" class="searchable">
    <table class="code"><thead><tr><th class="lineno" title="Line numbers">Line</th><th class="content"> </th></tr></thead><tbody><tr><th id="L1"><a href="#L1">1</a></th><td><span class="cp">#include &lt;libplugin/plugin.h&gt;</span></td></tr><tr><th id="L2"><a href="#L2">2</a></th><td><span class="cp">#include "psi4-dec.h"</span></td></tr><tr><th id="L3"><a href="#L3">3</a></th><td><span class="cp">#include &lt;libparallel/parallel.h&gt;</span></td></tr><tr><th id="L4"><a href="#L4">4</a></th><td><span class="cp">#include &lt;liboptions/liboptions.h&gt;</span></td></tr><tr><th id="L5"><a href="#L5">5</a></th><td><span class="cp">#include &lt;libmints/mints.h&gt;</span></td></tr><tr><th id="L6"><a href="#L6">6</a></th><td><span class="cp">#include &lt;libpsio/psio.h&gt;</span></td></tr><tr><th id="L7"><a href="#L7">7</a></th><td><span class="cp">#include &lt;libchkpt/chkpt.hpp&gt;</span></td></tr><tr><th id="L8"><a href="#L8">8</a></th><td><span class="cp">#include &lt;libqt/qt.h&gt;</span></td></tr><tr><th id="L9"><a href="#L9">9</a></th><td><span class="cp">#include "mad_mp2.h"</span></td></tr><tr><th id="L10"><a href="#L10">10</a></th><td><span class="cp"></span></td></tr><tr><th id="L11"><a href="#L11">11</a></th><td>INIT_PLUGIN</td></tr><tr><th id="L12"><a href="#L12">12</a></th><td></td></tr><tr><th id="L13"><a href="#L13">13</a></th><td><span class="k">namespace</span> psi<span class="p">{</span> </td></tr><tr><th id="L14"><a href="#L14">14</a></th><td></td></tr><tr><th id="L15"><a href="#L15">15</a></th><td><span class="k">extern</span> <span class="s">"C"</span> <span class="kt">int</span></td></tr><tr><th id="L16"><a href="#L16">16</a></th><td>read_options<span class="p">(</span>std<span class="o">::</span>string name<span class="p">,</span> Options <span class="o">&amp;</span>options<span class="p">)</span></td></tr><tr><th id="L17"><a href="#L17">17</a></th><td><span class="p">{</span></td></tr><tr><th id="L18"><a href="#L18">18</a></th><td>    <span class="k">if</span> <span class="p">(</span>name <span class="o">==</span> <span class="s">"MAD_MP2"</span><span class="o">||</span> options<span class="p">.</span>read_globals<span class="p">())</span> <span class="p">{</span></td></tr><tr><th id="L19"><a href="#L19">19</a></th><td>        <span class="cm">/*- The amount of information printed</span></td></tr><tr><th id="L20"><a href="#L20">20</a></th><td><span class="cm">            to the output file -*/</span></td></tr><tr><th id="L21"><a href="#L21">21</a></th><td>        options<span class="p">.</span>add_int<span class="p">(</span><span class="s">"PRINT"</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></td></tr><tr><th id="L22"><a href="#L22">22</a></th><td>        <span class="cm">/*- The amount of information printed</span></td></tr><tr><th id="L23"><a href="#L23">23</a></th><td><span class="cm">            to the output file -*/</span></td></tr><tr><th id="L24"><a href="#L24">24</a></th><td>        options<span class="p">.</span>add_int<span class="p">(</span><span class="s">"DEBUG"</span><span class="p">,</span> <span class="mi">1</span><span class="p">);</span></td></tr><tr><th id="L25"><a href="#L25">25</a></th><td>        <span class="cm">/*- The schwarz cutoff -*/</span></td></tr><tr><th id="L26"><a href="#L26">26</a></th><td>        options<span class="p">.</span>add_double<span class="p">(</span><span class="s">"SCHWARZ_CUTOFF"</span><span class="p">,</span> <span class="mf">0.0</span><span class="p">);</span></td></tr><tr><th id="L27"><a href="#L27">27</a></th><td>        <span class="cm">/*- The same-spin scale factor -*/</span></td></tr><tr><th id="L28"><a href="#L28">28</a></th><td>        options<span class="p">.</span>add_double<span class="p">(</span><span class="s">"SCALE_SS"</span><span class="p">,</span> <span class="mf">1.0</span><span class="o">/</span><span class="mf">3.0</span><span class="p">);</span></td></tr><tr><th id="L29"><a href="#L29">29</a></th><td>        <span class="cm">/*- The opposite-spin scale factor -*/</span></td></tr><tr><th id="L30"><a href="#L30">30</a></th><td>        options<span class="p">.</span>add_double<span class="p">(</span><span class="s">"SCALE_OS"</span><span class="p">,</span> <span class="mf">6.0</span><span class="o">/</span><span class="mf">5.0</span><span class="p">);</span></td></tr><tr><th id="L31"><a href="#L31">31</a></th><td>        <span class="cm">/*- Denominator algorithm for PT methods -*/</span></td></tr><tr><th id="L32"><a href="#L32">32</a></th><td>        options<span class="p">.</span>add_str<span class="p">(</span><span class="s">"DENOMINATOR_ALGORITHM"</span><span class="p">,</span> <span class="s">"LAPLACE"</span><span class="p">,</span> <span class="s">"LAPLACE CHOLESKY"</span><span class="p">);</span></td></tr><tr><th id="L33"><a href="#L33">33</a></th><td>        <span class="cm">/*- Maximum denominator error allowed (Max error norm in Delta tensor) -*/</span></td></tr><tr><th id="L34"><a href="#L34">34</a></th><td>        options<span class="p">.</span>add_double<span class="p">(</span><span class="s">"DENOMINATOR_DELTA"</span><span class="p">,</span> <span class="mf">1.0E-6</span><span class="p">);</span></td></tr><tr><th id="L35"><a href="#L35">35</a></th><td>        <span class="cm">/*- MP2 Algorithm -*/</span></td></tr><tr><th id="L36"><a href="#L36">36</a></th><td>        options<span class="p">.</span>add_str<span class="p">(</span><span class="s">"MP2_ALGORITHM"</span><span class="p">,</span> <span class="s">"DFMP2"</span><span class="p">,</span> <span class="s">"DFMP2 DFMP2J"</span><span class="p">);</span></td></tr><tr><th id="L37"><a href="#L37">37</a></th><td>    <span class="p">}</span></td></tr><tr><th id="L38"><a href="#L38">38</a></th><td>    <span class="k">return</span> <span class="kc">true</span><span class="p">;</span></td></tr><tr><th id="L39"><a href="#L39">39</a></th><td><span class="p">}</span></td></tr><tr><th id="L40"><a href="#L40">40</a></th><td></td></tr><tr><th id="L41"><a href="#L41">41</a></th><td><span class="k">extern</span> <span class="s">"C"</span> PsiReturnType</td></tr><tr><th id="L42"><a href="#L42">42</a></th><td>plugin_mad_mp2<span class="p">(</span>Options <span class="o">&amp;</span>options<span class="p">)</span></td></tr><tr><th id="L43"><a href="#L43">43</a></th><td><span class="p">{</span></td></tr><tr><th id="L44"><a href="#L44">44</a></th><td>    tstart<span class="p">();</span></td></tr><tr><th id="L45"><a href="#L45">45</a></th><td></td></tr><tr><th id="L46"><a href="#L46">46</a></th><td>    boost<span class="o">::</span>shared_ptr<span class="o">&lt;</span>PSIO<span class="o">&gt;</span> psio <span class="o">=</span> PSIO<span class="o">::</span>shared_object<span class="p">();</span></td></tr><tr><th id="L47"><a href="#L47">47</a></th><td></td></tr><tr><th id="L48"><a href="#L48">48</a></th><td>    boost<span class="o">::</span>shared_ptr<span class="o">&lt;</span>psi<span class="o">::</span>mad_mp2<span class="o">::</span>MAD_MP2<span class="o">&gt;</span> mp2<span class="p">(</span><span class="k">new</span> psi<span class="o">::</span>mad_mp2<span class="o">::</span>MAD_MP2<span class="p">(</span>options<span class="p">,</span> psio<span class="p">));</span></td></tr><tr><th id="L49"><a href="#L49">49</a></th><td>    </td></tr><tr><th id="L50"><a href="#L50">50</a></th><td>    mp2<span class="o">-&gt;</span>compute_energy<span class="p">();</span></td></tr><tr><th id="L51"><a href="#L51">51</a></th><td></td></tr><tr><th id="L52"><a href="#L52">52</a></th><td>    tstop<span class="p">();</span></td></tr><tr><th id="L53"><a href="#L53">53</a></th><td></td></tr><tr><th id="L54"><a href="#L54">54</a></th><td>    <span class="k">return</span> Success<span class="p">;</span></td></tr><tr><th id="L55"><a href="#L55">55</a></th><td><span class="p">}</span></td></tr><tr><th id="L56"><a href="#L56">56</a></th><td></td></tr><tr><th id="L57"><a href="#L57">57</a></th><td><span class="p">}</span> <span class="c1">// End Namespaces</span></td></tr></tbody></table>
      </div>
      <div id="help">
        <strong>Note:</strong> See <a href="/trac/wiki/TracBrowser">TracBrowser</a>
        for help on using the browser.
      </div>
      <div id="anydiff">
        <form action="/trac/diff" method="get">
          <div class="buttons">
            <input type="hidden" name="new_path" value="/trunk/tests/plugin_mad_mp2/main.cc" />
            <input type="hidden" name="old_path" value="/trunk/tests/plugin_mad_mp2/main.cc" />
            <input type="hidden" name="new_rev" />
            <input type="hidden" name="old_rev" />
            <input type="submit" value="View changes..." title="Select paths and revs for Diff" />
          </div>
        </form>
      </div>
    </div>
    <div id="altlinks">
      <h3>Download in other formats:</h3>
      <ul>
        <li class="first">
          <a rel="nofollow" href="/trac/browser/trunk/tests/plugin_mad_mp2/main.cc?format=txt">Plain Text</a>
        </li><li class="last">
          <a rel="nofollow" href="/trac/export/1983/trunk/tests/plugin_mad_mp2/main.cc">Original Format</a>
        </li>
      </ul>
    </div>
    </div>
    <div id="footer" lang="en" xml:lang="en"><hr />
      <a id="tracpowered" href="http://trac.edgewall.org/"><img src="/trac/chrome/common/trac_logo_mini.png" height="30" width="107" alt="Trac Powered" /></a>
      <p class="left">
        Powered by <a href="/trac/about"><strong>Trac 0.11.7</strong></a><br />
        By <a href="http://www.edgewall.org/">Edgewall Software</a>.
      </p>
      <p class="right">www.psicode.org - The home of Psi4 on the internet.</p>
    </div>
  </body>
</html>